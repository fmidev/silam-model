"""

Two functions for creating sets of text files based on a template
and replacement rules. Example:

Assume file 'foo.in' contains strings $DATE and $NUM. We create 10
files corresponding to first days of years 2300-2309. 

import datetime
from file_generator import generate
times = [datetime.datetime(i, 1, 1) for i in range(2300,2310)]
substitutions = {'$DATE' : [t.strftime('%F') for t in times], 
                 '$NUM'  : range(10) }
outfiles = ['foo_%i.out' % i for i in range(10)]
generate('foo.in', outfiles, substitutions)

---
Also the following function expand could have been used as 
... = {'$DATE' : expand('%F', times), ...

"""


def generate(file_template, files_out, substitutions):
    # Input arguments: 
    # 
    # file_template : path to the template file
    # files_out : sequence of output file names 
    # substitutions : dictionary {key : replacements}, where key is the replaced
    # string and replacements is a sequence. The replacements are converted into 
    # strings, if possible.
    
    # The replacement lists are required to have the same length as files_out.

    values = substitutions.values()

    if isinstance(files_out, str):
        files_out = (files_out,)
        nfiles = 1
        if not all([len(x) == len(values[0]) for x in values]):
            raise ValueError('Substitution vectors must have same length')
    else:
        nfiles = len(files_out)
        if not all([len(x) == nfiles for x in values]):
            raise ValueError('Substitution vectors must have same length')
    
    
    

    f_in = open(file_template, 'r')
    
    if nfiles == 1:
        f_out = open(files_out[0], 'w')

    for i in range(len(values[0])):
        if nfiles > 1:
            f_out = open(files_out[i], 'w')
        
        for line_in in f_in:
            line = line_in
            for (key, val) in substitutions.items():
                if key in line_in:
                    line = line.replace(key, str(val[i]))

            f_out.write(line)
    
        if nfiles > 1:
            f_out.close()
        f_in.seek(0,0)


def expand(template, times=None, nums=None):
    if times and nums:
        if not len(times) == len(nums):
            raise ValueError('Times and nums must have same length')
        return [t.strftime(template).replace('#', str(n)) for (t,n) in zip(times, nums)]
    
    if times:
        return [t.strftime(template) for t in times]
    
    if nums:
        return [template.replace('#', str(n)) for n in nums]

        
            
    


                 
