'''
Created on 18 Jul 2013

@author: vimaier

GetLLM.test.vimaier_utils.compare_utils holds functions for comparing TfsFiles.

'''
import string


import metaclass

#===================================================================================================
# helper-functions
#===================================================================================================
def parse_descriptor_line(line):
    '''
    Parses a TfsFile decriptor line('@ XYZ %s "A_string"').
    Returns the single parts in a list(['@', 'XYZ', '%s', 'A_string']
    
    :Parameters:
        'line': string
            TfsFile descriptor line
    :Return: list
        List with single parts
    '''
    start_i, end_i = __get_index_of_name(line)
    descriptor_name = line[start_i:end_i]
    
    start_i, end_i = __get_index_of_type(line, end_i)
    descriptor_type = line[start_i:end_i]

    start_i, end_i = __get_index_of_value(line, end_i)
    descriptor_value = line[start_i:end_i]
    
    return ["@", descriptor_name, descriptor_type, descriptor_value]

def __get_index_of_name(line):
    """ Returns for example string '@ xyz 5s....' (2,5) """
    start_i = 0
    end_i = 0
    
    # Start index
    while True:
        if is_whitespace(line[start_i]) or '@' == line[start_i]:
            start_i += 1
            continue
        else:
            break
    
    # End index
    for i in xrange(start_i+1, len(line)):
        if is_whitespace(line[i]):
            end_i = i
            break
    
    return start_i, end_i


def __get_index_of_type(line, start_i):
    #Find index of '%'
    start_i = str.find(line, '%', start_i)
    
    # End index
    for i in xrange(start_i+1, len(line)):
        if is_whitespace(line[i]):
            end_i = i
            break
    
    return start_i, end_i

def __get_index_of_value(line, start_i):
    # Find first character
    for i in xrange(start_i, len(line)):
        if not is_whitespace(line[i]):
            start_i = i
            break
        
    if '"' == line[start_i]:
        # Find next '"'
        end_i = str.find(line, '"', start_i+1)
        # We don't want the character '"' in the string:
        start_i += 1
    else:
        # Find next whitespace
        end_i = len(line)
        for i in xrange(start_i+1, len(line)):
            if is_whitespace(line[i]):
                end_i = i
                break
    
    return start_i, end_i
    
def is_whitespace(single_char):
    for whitespace in string.whitespace:
        if single_char == whitespace:
            return True
    
    return False



def compare_tfs_files(name_valid, name_to_check):
    """ Compares both files. Whitespace does not matter.
        Returns an error message or in success an empty string.
    """
    if name_valid.endswith(".gitignore"):
        return ""
    
    if __file_not_valid(name_valid):
        return ""
    
    file_valid = open(name_valid)
    file_to_check = open(name_to_check)
    
    valid_lines = file_valid.readlines()
    to_check_lines = file_to_check.readlines()
    
    i_to_check = 0
    for i_valid in xrange(len(valid_lines)):
        valid_line = valid_lines[i_valid]
        
        # Exclude descriptors GetLLMVersion, MAD_FILE, FILES, DATE and Command from comparison
        if valid_line.startswith("@") and "GetLLMVersion" in valid_line or "MAD_FILE" in valid_line or "FILE" in valid_line:
            continue
        while to_check_lines[i_to_check].startswith("@") and "GetLLMVersion" in to_check_lines[i_to_check] or "MAD_FILE" in to_check_lines[i_to_check] or "FILE" in to_check_lines[i_to_check] or "DATE" in to_check_lines[i_to_check] or "Command" in to_check_lines[i_to_check]:
            i_to_check += 1
        check_line = to_check_lines[i_to_check]
            
        if valid_line.startswith("@"):
            split_valid = parse_descriptor_line(valid_line)
        else:            
            split_valid = valid_line.split()
        
        if check_line.startswith("@"):
            split_to_check = parse_descriptor_line(check_line)
        else:
            split_to_check = check_line.split()

        if len(split_valid) != len(split_to_check):
            return "Column numbers not equal:\n"+valid_line+check_line
        
        for i in xrange(len(split_valid)):
            if split_valid[i] != split_to_check[i]:
                err_msg = "Entry in column number["+str(i)+"]not equal:\n"+valid_line+check_line
                err_msg += str(split_valid[i]) +" != "+ str(split_to_check[i])
                return err_msg
        i_to_check += 1
    
    return ""
        

def __file_not_valid(name_twiss):
    ''' Checks if the given twiss file is empty, thus not valid. '''
    try:
        tw = metaclass.twiss(name_twiss)
        if 0 == len(getattr(tw, "NAME", [])):
            return True # empty files are not longer produced by modified GetLLM
        return False
    except ValueError:
        # Probably empty file
        return True