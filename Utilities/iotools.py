'''
Created on 1 Jul 2013

@author: vimaier

@version: 1.0.0

Utilities.iotools.py holds helper functions for input/output issues. This module is not intended to 
be executed.

Feel free to use and extend this module.

'''

import os
import sys
import shutil




def delete_content_of_dir(path_to_dir):
    ''' 
    Deletes all folders, files and symbolic links in given directory.
    :parameters:
        path_to_dir:string
    '''
    if not os.path.isdir(path_to_dir):
        return
    
    for item in os.listdir(path_to_dir):
        item_path = os.path.join(path_to_dir, item)
        delete_item(item_path)


def delete_item(path_to_item):
    ''' Deletes the item given by path_to_item. It distinguishes between a file, a directory and a
    symbolic link.
    '''
    try:
        if os.path.isfile(path_to_item):
            os.unlink(path_to_item)
        elif os.path.isdir(path_to_item):
            shutil.rmtree(path_to_item)
        elif os.path.islink(path_to_item):
            os.unlink(path_to_item)
    except IOError:
        print >> sys.stderr, "Could not delete item because of IOError.", path_to_item
        
        
def copy_content_of_dir(src_dir, dst_dir):
    ''' Copies all files and directories from src_dir to dst_dir. '''
    if not os.path.isdir(src_dir):
        return
    
    create_dirs(dst_dir)
    
    for item in os.listdir(src_dir):
        src_item = os.path.join(src_dir, item)
        dst_item = os.path.join(dst_dir, item)
        copy_item(src_item, dst_item)
                

def create_dirs(path_to_dir):
    ''' Creates all dirs to path_to_dir if not exists. '''
    if not os.path.exists(path_to_dir):
        os.makedirs(path_to_dir)
    
    
def copy_item(src_item, dest):
    ''' Copies a file or a directory to dest. dest may be a directory.
    If src_item is a directory then all containing files and dirs will be copied into dest. '''
    try:
        if os.path.isfile(src_item):
            shutil.copy2(src_item, dest)
        elif os.path.isdir(src_item):
            copy_content_of_dir(src_item, dest)
    except IOError:
        print >> sys.stderr, "Could not copy item because of IOError.", src_item
        
        
def deleteFilesWithoutGitignore(pathToDirectory):
    """
    Deletes all files in the given pathToDirectory except of the file with the name '.gitignore'
    
    :Return: boolean
          True if the directory exists and the files are deleted otherwise False. 
    """
    if not os.path.exists(pathToDirectory):
        return False
    
    filenames_list = os.listdir(pathToDirectory)
    
    for filename in filenames_list:
        if ".gitignore" != filename:
            os.remove( os.path.join(pathToDirectory,filename) )
    
    return True

def existsDirectory(path_to_dir):
    return os.path.isdir(path_to_dir)

def notExistsDirectory(path_to_dir):
    return not existsDirectory(path_to_dir)

def get_absolute_path_to_betabeat_root():
    print os.path.dirname(os.path.abspath(__file__))
    print os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.pardir)
    return os.path.abspath(
                    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.pardir)
                    )
    
    