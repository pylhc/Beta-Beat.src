r"""
.. module: Utilities.convert2json

Created on ??

Loads the dictionaries which are the return values of the AllLists.py python script
functions and dumps it as a json file.

.. moduleauthor:: Yngve
"""

import json
import os


def convert(fin, fout):
    print "Converting", fin, "-->", fout
    result_dict = {}
    module_name = fin[:-3]
    list_of_module_attributes = []
    exec('import '+module_name)
    exec('list_of_module_attributes = dir('+module_name+')')
    for key in list_of_module_attributes:
        if '_' not in key:
            print '  ', key
            exec('result_dict[key] = '+module_name+'.'+key+'()')
    file(fout,'w').write(json.dumps(result_dict, sort_keys=True, indent=2))

if __name__ == "__main__":
    for f in os.listdir('.'):
        if f[:8] == 'AllLists' and f[-3:] == '.py':
            fout = f[:-2] + 'json'
            fin = f
            convert(fin, fout)


