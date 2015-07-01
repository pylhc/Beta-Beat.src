import __init__  # @UnusedImport
import sys
import os
import re
from Utilities import iotools

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class MadxTemplates():
    """
    This class is meant to load all the templates in the "templates" directory,
    together with those in the input "template_dirs" directory. It will generate
    a method in the class for each template, with the placeholders as attributes.
    When those functions are called, "madx_wrapper" will run replacing the attributes
    into the template.
    """

    def __init__(self, template_dirs=[], output_file=None, log_file=None, verbose=False):
        self.__output_file = output_file
        self.__log_file = log_file
        self._verbose = verbose
        template_dirs.append(os.path.join(CURRENT_PATH, "templates"))
        if self._verbose:
            print "Using template directories: " + str(template_dirs)
            print "Generating methods..."
        for template_dir in template_dirs:
            if os.path.isdir(template_dir):
                for template_name in os.listdir(template_dir):
                    template_path = os.path.join(template_dir, template_name)
                    if os.path.isfile(template_path):
                        placeholder_list = self.__parse_template(template_path)
                        if len(placeholder_list) != 0:
                            self.__assing_new_method(template_name, template_path, placeholder_list)
                        elif self._verbose:
                            print "WARNING: Ignoring non-template: " + os.path.abspath(template_path) + ". Cannot find placeholders."
                    elif self._verbose:
                        print "WARNING: Ignoring non-file: " + os.path.abspath(template_path)
            elif self._verbose:
                print "WARNING: Ignoring non-directory: " + os.path.abspath(template_dir)
        if self._verbose:
            print "All done."

    def set_log_path(self, log_path):
        self.log_file = log_path

    def set_output_path(self, output_path):
        self.output_file = output_path

    def __parse_template(self, template_path):
        template_content = iotools.read_all_lines_in_textfile(template_path)
        placeholder_list = self.__get_placeholders_from_string(template_content)
        return placeholder_list

    def __get_placeholders_from_string(self, string):
        placeholder_list = re.findall("%\((.*?)\)", string)
        return self.__remove_duplicated(placeholder_list)

    def __remove_duplicated(self, item_list):
        new_list = []
        for item in item_list:
            if item not in new_list:
                new_list.append(item)
        return new_list

    def __assing_new_method(self, template_name, template_path, placeholder_list):
        function_name = template_name.replace(".", "_")
        new_function_code = "def " + function_name + "("
        placeholder_replacement_dict = "{"
        for placeholder in placeholder_list:
            new_function_code += placeholder + ", "
            placeholder_replacement_dict += "\"" + placeholder + "\"" + ": " + placeholder + ", "
        placeholder_replacement_dict = placeholder_replacement_dict[:-2] + "}"
        new_function_code = new_function_code[:-2] + "):\n"
        if self._verbose:
            print new_function_code[:-2]
        output_file = None if self.__output_file is None else "\"" + self.__output_file + "\""
        log_file = None if self.__log_file is None else "\"" + self.__log_file + "\""
        new_function_code += "    from Utilities import iotools\n"
        new_function_code += "    from madx import madx_wrapper\n"
        new_function_code += "    template_content = iotools.read_all_lines_in_textfile(r\"" + template_path + "\")" + "\n"
        new_function_code += "    resolved_template = template_content % " + placeholder_replacement_dict + "\n"
        new_function_code += "    madx_wrapper.resolve_and_run_string(resolved_template, output_file=" + str(output_file) + ", log_file=" + str(log_file) + ")\n"
        exec new_function_code in self.__dict__


if __name__ == '__main__':
    MadxTemplates(verbose=True)
    print >> sys.stderr, "This module is only meant to be imported."
    sys.exit(-1)
