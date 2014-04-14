import sys
import json


ALL_LISTS_BEAM1_PATH = '/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB/fullresponse/LHCB1/AllLists.json'
ALL_LISTS_BEAM2_PATH = '/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB/fullresponse/LHCB2/AllLists.json'


def main():
    variables_beam1 = json.load(file(ALL_LISTS_BEAM1_PATH, 'r'))['getListsByIR'][1]
    variables_common, variables_beam2 = json.load(file(ALL_LISTS_BEAM2_PATH, 'r'))['getListsByIR']

    ip_int = 1
    if len(sys.argv) > 1:
        ip_int = int(sys.argv[1])
    if ip_int > 8 or ip_int < 1:
        ip_int = 1
    ip_string = str(ip_int)

    apply_correction_file = open("applycorrection.seqx", 'w')
    variables_common_file = open("variablesc.seqx", 'w')
    variables_beam1_file = open("variablesb1.seqx", 'w')
    variables_beam2_file = open("variablesb2.seqx", 'w')
    variables_s_file = open("svariables.seqx", 'w')
    variables_d_file = open("dvariables.seqx", 'w')

    param_change_generator_file = open("genchangpars.seqx", 'w')
    param_change_generator_file.write('select,flag=save, clear;')

    variables = variables_beam1[ip_string]
    print '\nBeam 1\n'
    param_change_generator_file.write('!B1\n')
    _vars_to_files(apply_correction_file, variables_beam1_file, variables_s_file, variables_d_file, param_change_generator_file, variables)

    variables = variables_beam2[ip_string]
    print '\nBeam 2\n'
    param_change_generator_file.write('\n!B2\n')
    _vars_to_files(apply_correction_file, variables_beam2_file, variables_s_file, variables_d_file, param_change_generator_file, variables)

    variables = variables_common[ip_string]
    print '\nBeam 1 and Beam 2\n'
    param_change_generator_file.write('\n!B1 and B2\n')
    _vars_to_files(apply_correction_file, variables_common_file, variables_s_file, variables_d_file, param_change_generator_file, variables)

    variables_common_file.close()
    variables_beam1_file.close()
    variables_beam2_file.close()
    variables_s_file.close()
    variables_d_file.close()

    param_change_generator_file.write('\n save, file=\"changeparameters.madx\";\n')
    param_change_generator_file.close()


def _vars_to_files(apply_correction_file, variables_file, variables_s_file, variables_d_file, param_change_generator_file, variables):
    for variable in variables:
        print variable
        variables_file.write('   vary, name=d' + variable + ', step:=1e-4;\n')
        variables_s_file.write(' ' + variable + '_0 = ' + variable + ';\n')
        variables_d_file.write(' ' + variable + ' := ' + variable + '_0 + d' + variable + ';\n')
        param_change_generator_file.write('select,flag=save,pattern=\"d' + variable + '\";\n')
        apply_correction_file.write(variable + ' = ' + variable + '_0 + d' + variable + ';\n')


if __name__ == "__main__":
    main()
