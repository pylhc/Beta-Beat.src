import sys
import subprocess
import logging
''' Wrapper for jws to execute java programs distributed as jnlp
    
    author: Piotr Skowronski
    
    Usage ::
        python jws_wrapper.py pathToJws JwsOpton1 JwsOpton2 ... JnlpUrl JnlpOption1 JnlpOption2 ... 
      
    This script mangles the parameters and executes via system call jws like this
    pathToJws JwsOpton1 JwsOpton2 ... JnlpUrl?arg1=JnlpOption1&arg2=JnlpOption2& ...

'''

#logging.basicConfig()
LOGGER = logging.getLogger("JwsWrapper")

    
class JwsWrapper:
    def __init__(self):
        if sys.flags.debug: # true of -d passed to python
            LOGGER.setLevel(logging.DEBUG)
            ch = logging.StreamHandler(sys.stdout)
            ch.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            ch.setFormatter(formatter)
            LOGGER.addHandler(ch)
            
        self.jws_path = ""
        self.jws_pars = []
        self.jnlp_url = ""
        self.jnlp_pars = "" # arg1=-e&arg2=-time&arg3=...
        
    def mangleargs(self,argv):
        ''' Decodes arguments.  
            Argument 1 is jws path
            The argument ending with .jnlp is the url of the jnlp to execute.
            All arguments before are options of jws
            All arguments after are options of jnlp
        '''
        
        LOGGER.debug("This is the name of the script: %s ", argv[0])
        LOGGER.debug("Number of arguments: %d", len(argv))
        LOGGER.debug("The arguments are: %s" , str(argv))
        
        if len(argv) < 3:
            raise SyntaxError('Needs 2 arguments: path to jws and a jnlp url')
        self.jws_path = argv[1]
        
        idxjnlp=0
        for i in range(2,len(argv)):
            a = argv[i]
            if a.endswith(".jnlp"):
                idxjnlp = i
                break
            else:
                self.jws_pars.append(a)
        
        if idxjnlp < 1:
            raise SyntaxError('Did not find a jnlp in list of arguments')
        
        self.jnlp_url = argv[idxjnlp]
        
        # arg1=-e&arg2=-time&arg3=...
        for i in range(idxjnlp+1,len(argv)):
            
            argno = i - idxjnlp
            LOGGER.debug("jnlp arg no %d %s ", argno, argv[i])

            argstr = "arg{}={}".format(argno,argv[i])
            if argno > 1:
                argstr = "&"+argstr
            
            argstr = argstr.replace(" ", "_")
            self.jnlp_pars += argstr
            
        
        return 0
    def calljws(self):
        #self.jws_path = "ls"
        #self.jws_pars.append("-l")
        #self.jws_pars.append("/home")
        #self.jnlp_url = "blabla.jnlp"
        #self.jnlp_pars = "&par1=-m&par2=-e"

        LOGGER.debug( "JWS path: %s", self.jws_path)
        LOGGER.debug( "JWS options: %s", self.jws_pars)
        LOGGER.debug( "JNLP URL: %s" , self.jnlp_url)
        LOGGER.debug( "JNLP args: %s" , self.jnlp_pars)

        
        self.jws_pars.append(self.jnlp_url +"?" +self.jnlp_pars)
        subprocess.call([self.jws_path] + self.jws_pars)
        
    def runjws(self, argv):
        
        self.mangleargs(argv)
        retval =  self.calljws()
        
        return retval
        
               
def runjws(argv):
    jwsw = JwsWrapper()
    retval = jwsw.runjws(argv)
    return retval
    

    
if __name__ == "__main__":
    return runjws(sys.argv)
