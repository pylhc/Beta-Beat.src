'''
Created on 22 May 2013

@author: vimaier

@version: 1.0.0

This module stores functions for input/output

Change history:

'''
import os


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
        
            
    
    