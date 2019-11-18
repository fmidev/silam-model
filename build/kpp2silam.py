"""
This script creates a SILAM interface to a KPP generated chemistry solver. This includes
substituting a number of flags into the _template files included in silam v5.3
source. Special steps are taken to make the KPP module variables OpenMP threadprivate.

After running this script, the chemistry modules should be almost usable - the only manual
changes should be in the top-level chem_dep file. This is why it is not overwritten by
default. Obviously, the new chemistry scheme needs still to be introduced in chemistry
server, manually.
"""

import os, argparse
from os import path
import file_generator as fgen
import shutil
import re

parser = argparse.ArgumentParser('kpp2silam')
parser.add_argument('--kpp-dir', '-i', help='path to directory with the .kpp file', default='.')
parser.add_argument('--silam-dir', '-o', help='SILAM source directory', default='.')
parser.add_argument('kpp_root', help='Name of the scheme (in both kpp and SILAM)')
parser.add_argument('--overwrite', '-f', action='store_true', default=False, 
                    help='Overwrite existing chem_dep files')
parser.add_argument('--transf-id', default=8888,
                    help='Substitute transformation number, chem_dep module', type=int)
parser.add_argument('--use-module', '-u', help="Include use statement for this module in interface",
                    action='append')
parser.add_argument('--photoarr', action='store_true', default=False,
                    help='Include the photoarr argument to set_rates subroutine')
parser.add_argument('--SOA', action='store_true', default=False,
                    help='Include the secondary organic aerosols')
args = parser.parse_args()

kpp_root = args.kpp_root
kpp_dir = args.kpp_dir
kpp_file = path.join(kpp_dir, '%s.kpp' % kpp_root)
silam_src_dir = args.silam_dir
interface_template = path.join(silam_src_dir, 'kpp_interface')
chem_dep_template = path.join(silam_src_dir, 'chem_dep_template')
kpp_global_template = path.join(silam_src_dir, 'kpp_global_template')

if not 'KPP_HOME' in os.environ or os.environ['KPP_HOME'].isspace():
    os.environ['KPP_HOME'] = path.join(kpp_dir, '..')
print os.environ['KPP_HOME']

substitutions = {'$KPP_ROOT' : [kpp_root], '$TRANSF_ID' : [args.transf_id]}

def get_filename(which):
    return path.join(kpp_dir, '%s_%s.f90' % (kpp_root, which))

#modules_copied = 'Integrator Jacobian JacobianSP LinearAlgebra Parameters Precision'.split()
modules_copied = 'Integrator Jacobian JacobianSP LinearAlgebra Parameters'.split()
files_copied = [get_filename(module) for module in modules_copied]

def parse_kpp_parameters():
    filename_par = get_filename('Parameters')
    subst_fixed = {}
    subst_var = {}
    reader = open(filename_par, 'r')
    for line in reader:
        if 'Index declaration' in line:
            break
    reader.next()
    reader.next()
    # variable species
    pattern_var = ':: +ind_(.+) = (.+)'
    for line in reader:
        if line.isspace():
            break
        subst, ind = re.search(pattern_var, line).groups()
        ind = int(ind)
        subst_var[ind] = subst
    while not 'Index declaration for fixed species in FIX' in line:
        line = reader.next()

    pattern_fix = ':: +indf_(.+) = (.+)'
    reader.next()
    reader.next()
    for line in reader:
        if line.isspace():
            break
        subst, ind = re.search(pattern_fix, line).groups()
        ind = int(ind)
        subst_fixed[ind] = subst
    reader.close()
    subst_list_var = [subst_var[key] for key in sorted(subst_var.keys())]
    subst_list_fix = [subst_fixed[key] for key in sorted(subst_fixed.keys())]    
    return subst_list_var, subst_list_fix

def get_set_rates():
    filename = get_filename('Rates')
    reader = open(filename, 'r')
    line = ''
    lines = []
    indent = '    '
    while not 'Update_RCONST' in line:
        line = reader.next()
    for line in reader:
        if line.startswith('END'):
            break
        if re.search('RCONST\(.+\) =', line):
            lines.append(indent + line.strip())
    #while not 'Update_PHOTO' in line:
    #    line = reader.next()
    #for line in reader:
    #    if line.startswith('END'):
    #        break
    #    if re.search('RCONST\(.+\) =', line):
    #        lines.append(indent + line.strip())
    reader.close()
    return lines

def get_set_const_rates():
    indent = '    '
    filename = get_filename('Initialize')
    reader = open(filename, 'r')
    lines = []
    for line in reader:
        if re.search('RCONST\(.+\) =', line):
            lines.append(indent + line.strip())
    reader.close()
    return lines

def get_subst_names(subst_name_list):
    num_cols = 6
    ind_row = 0
    ind_col = 0
    lines = []
    iter_subst = iter(subst_name_list)
    max_subst_name_len = max(len(s) for s in subst_name_list)
    # fortran wants the name list to consist of elements with same length
    fmt = """'%%-%is'""" % max_subst_name_len
    out_of_values = False
    while True:
        lines.append([])
        for count_col in range(num_cols):
            try:
                elem = fmt % iter_subst.next()
                lines[-1].append(elem)
            except StopIteration:
                out_of_values = True
                break
        if out_of_values:
            break
    strlines = [','.join(line) for line in lines if len(line) > 0]
    multiline = ', &\n& '.join(strlines)
    return multiline
        
            

    
    # out_of_values = False
    # while True:
    #     if lines:
    #         lines[-1] += ', &'
    #     line = []
    #     try:
    #         for ind_col in range(num_cols):
    #             subst = iter_subst.next()
    #             subst = fmt % subst
    #             line.append("'%s'" % subst)
    #     except StopIteration:
    #         out_of_values = True
    #     finally:
    #         lines.append('& ' + ','.join(line))
    #         ind_row += 1
    #     if out_of_values:
    #         break
    # return '\n'.join(lines)
        
def cp_insrt_threadprivate(path_in, path_out, variables_thr_priv, var_type='real'):
    read_in = open(path_in, 'r')
    write_out = open(path_out, 'w')
    pattern_declr = '%s.*::' % var_type
    pattern_varname = '(\w+)(:?\(\w+\))?'
    pattern_var_declr = '(%s.*'
    header_done = False
    for line in read_in:
        line_lower_lstrip = line.lower().lstrip()
        if line_lower_lstrip.startswith('contains'):
            header_done = line.lstrip().lower().startswith('contains')
        if line_lower_lstrip.startswith('use'):
            write_out.write('!$use omp_lib\n')
        #if header_done:
        #    write_out.write(line)
        #    continue
        # while still in header, scan the module variables
        varnames = []
        if not header_done and re.match(pattern_declr, line_lower_lstrip):
            declr, var_list = line.split('::')
            var_declrs = var_list.split(',')
            for var in var_declrs:
                varname = re.match(pattern_varname, var.replace(' ', '')).groups()[0]
                if varname in variables_thr_priv:
                    print 'Flagged threadprivate: %s in %s' % (varname, path_in)
                    varnames.append(varname)
        if varnames:
            write_out.write('%s, save :: %s\n' % (declr, var_list))
            write_out.write('!$OMP THREADPRIVATE(%s)\n' % ','.join(varnames))
        else:
            write_out.write(line)
    read_in.close()
    write_out.close()

def cp_change_uses(path_in, path_out):
    read_in = open(path_in, 'r')
    write_out = open(path_out, 'w')
    for line in read_in:
        if line.lower().lstrip().startswith('use'):
            line_out = line.replace('%s_adj' % kpp_root, kpp_root)
        else:
            line_out = line
        write_out.write(line_out)
    read_in.close()
    write_out.close()

def add_dp_def(path_in, path_out):
    read_in = open(path_in, 'r')
    write_out = open(path_out, 'w')
    for line in read_in:
        line_lower_lstrip = line.lower().lstrip()
        if line_lower_lstrip.startswith('integer, parameter :: sp = selected_real_kind(6,30)\n'):
            write_out.write('#ifdef DOUBLE_PRECISION\n')
            write_out.write('  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(14,300)\n')
            write_out.write('#else\n')
            write_out.write(line)
            write_out.write('#endif\n')
        else:
            write_out.write(line)
    read_in.close()
    write_out.close()
    
class KPPError(Exception):
    def __init__(self):
        Exception.__init__(self, 'KPP failed')


command = path.expandvars('cd %s && $KPP_HOME/bin/kpp %s' % (kpp_dir, path.basename(kpp_file)))
if os.system(command):
    raise KPPError()

subst_var, subst_fixed = parse_kpp_parameters()
fixed_args = ['cnc_%s' % subst for subst in subst_fixed]
lines_set_fixed = ['fix(%i) = %s' % (ind+1, arg) for ind, arg in enumerate(fixed_args)]
set_fixed = '\n'.join(lines_set_fixed)
substitutions['$SET_FIXED'] = [set_fixed]
substitutions['$FIXED_ARGS'] = [','.join(fixed_args)]

subst_names = ','.join(subst_var)

set_rates = '\n'.join(get_set_rates())
substitutions['$SET_RATES'] = [set_rates]

set_const_rates = '\n'.join(get_set_const_rates())
substitutions['$SET_CONST_RATES'] = [set_const_rates]

substitutions['$SUBST_NAMES'] = [get_subst_names(subst_var)]



kpp_file_adj = kpp_file.replace('.kpp', '_adj.kpp')
have_adjoint = path.exists(kpp_file_adj)
if have_adjoint:
    command = path.expandvars('cd %s && $KPP_HOME/bin/kpp %s' % (kpp_dir, path.basename(kpp_file_adj)))
    if os.system(command):
        raise KPPError()
    substitutions['%s_adj'%kpp_root] = [kpp_root]
    files_copied.append(get_filename('Hessian'))
    files_copied.append(get_filename('HessianSP'))
    #files_copied.append(get_filename('adj_Integrator'))
    substitutions['$COMMENT_IF_NO_ADJ'] = ['']
    substitutions['$COMMENT_IF_HAVE_ADJ'] = ['!']
    # make the adjoint integrator refer to the basic modules of the forward model.
    adj_integrator_out = path.join(silam_src_dir, '%s_adj_Integrator.f90'%kpp_root)
    cp_change_uses(get_filename('adj_Integrator'), adj_integrator_out)
else:
    substitutions['$COMMENT_IF_NO_ADJ'] = ['!']
    substitutions['$COMMENT_IF_HAVE_ADJ'] = ['']

if args.photoarr:
    substitutions['$COMMENT_IF_NO_PHOTOARR'] = ['']
    substitutions['$COMMENT_IF_PHOTOARR'] = ['!']
else:
    substitutions['$COMMENT_IF_NO_PHOTOARR'] = ['!']
    substitutions['$COMMENT_IF_PHOTOARR'] = ['']
    
if args.SOA:
    substitutions['$COMMENT_IF_NO_SOA'] = ['']
    substitutions['$COMMENT_IF_SOA'] = ['!']
else:
    substitutions['$COMMENT_IF_NO_SOA'] = ['!']
    substitutions['$COMMENT_IF_SOA'] = ['']
    
if args.use_module:
    substitutions['$EXTRA_USE_LINES'] = ['\n'.join(['use %s' % module for module in args.use_module])]
else:
    substitutions['$EXTRA_USE_LINES'] = ['']
    
interface_out = path.join(silam_src_dir, '%s_interface.mod.f90' % kpp_root)
fgen.generate(interface_template, interface_out, substitutions)

chem_dep_out = path.join(silam_src_dir, 'chem_dep_%s.mod.f90' % kpp_root)
exists = path.exists(chem_dep_out)
if not exists or args.overwrite:
    if exists:
        shutil.move(chem_dep_out, chem_dep_out + '.bak')
    fgen.generate(chem_dep_template, chem_dep_out, substitutions)

kpp_function_out = path.join(silam_src_dir, '%s_Function.f90' % kpp_root)
cp_insrt_threadprivate(get_filename('Function'), kpp_function_out, ['A'])

kpp_global_out = path.join(silam_src_dir, '%s_Global.f90' % kpp_root)
fgen.generate(kpp_global_template, kpp_global_out, substitutions)

kpp_precision_out = path.join(silam_src_dir, '%s_Precision.f90' % kpp_root)
add_dp_def(get_filename('Precision'), kpp_precision_out)

for filename in files_copied:
    shutil.copy(filename, silam_src_dir)

