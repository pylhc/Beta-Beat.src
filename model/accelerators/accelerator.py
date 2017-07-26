class AccExcitationMode(object):
    FREE, ACD, ADT = range(3)


class Accelerator(object):
    """
    Abstract class to serve as an interface to implement the
    rest of the accelerators.
    """

    # Class methods ###########################################

    @classmethod
    def get_class(cls, *args, **kwargs):
        """
        This method should return the accelerator class defined
        in the arguments.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def get_class_from_args(cls, args):
        """
        This method should return the accelerator class defined
        by in the passed command line like arguments.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def init_from_args(cls, args):
        """
        Instances an accelerator filling the attributes using the
        command line like arguments passed in args.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def get_variables(cls, frm=None, to=None, classes=None):
        """
        Gets the variables with elements in the given range and the given
        classes. None means everything.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def get_bpms_around_acd(available_bpm_list, model_twiss, exciter="acd"):
        """
        It will return the two best BMPs to use for exciter specified
        in compensation.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def get_arc_bpms_mask(list_of_elements):
        """
        It will return a mask to filter with the elements of @list_of_elements
        so that only arc BPMs remains.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def get_correctors_variables(cls, frm=None, to=None, classes=None):
        """
        Returns the set of corrector variables between frm and to, with classes
        in classes. None means select all.
        """
        Accelerator.__raise_definition_error()

    # Instance methods ########################################

    def verify_object(self):
        """
        Verifies that this instance of an accelerator is properly
        instantiated.
        """
        Accelerator.__raise_definition_error()

    @staticmethod
    def __raise_definition_error():
        raise AcceleratorDefinitionError(
            """A function that should have been overwritten
            has been called, check stack trace."""
        )
        
    # For GetLLM -------------------------------------------
    
    @classmethod
    def get_beam_direction(self):
        Accelerator.__raise_definition_error()
    
    def get_exciter_bpm(self, plane, distance):
        """
        Returns the BPM next to the exciter.
        The accelerator instance knows already which excitation method is used.
        distance: 1=nearest bpm 2=next to nearest bpm
        """
        Accelerator.__raise_definition_error()
    
    def get_exciter_name(self, plane):
        """
        Returns the name of the exciter.
        """
        Accelerator.__raise_definition_error()

    def get_model_tfs(self):  # instance method because it has to access the instance's model
        """
        Returns the model tfs file.
        """
        Accelerator.__raise_definition_error()

    def get_driven_tfs(self):
        """
        Returns the driven model tfs file.
        """
        Accelerator.__raise_definition_error()

    def get_best_knowledge_model_tfs(self):
        """
        Returns the best knowledge model tfs file.
        """
        Accelerator.__raise_definition_error()

    def get_elements_tfs(self):
        """
        Returns theelements tfs file.
        """
        Accelerator.__raise_definition_error()
        
    def get_s_first_BPM(self):
        """
        Returns the position of the first BPM in turn by turn acquisition.
        """

    ###########################################################


class Variable(object):
    """
    Generic corrector variable class that holds name, position (s) and
    physical elements it affectes. This variables should be logical variables
    that should have and effect in the model if modified.
    """
    def __init__(self, name, elements, classes):
        self.name = name
        self.elements = elements
        self.classes = classes


class Element(object):
    """
    Generic corrector element class that holds name and position (s)
    of the corrector. This element should represent a physical element of the
    accelerator.
    """
    def __init__(self, name, s):
        self.name = name
        self.s = s


class AcceleratorDefinitionError(Exception):
    """
    Raised when an accelerator instance is wrongly used, for
    example by calling a method that should have been overwritten.
    """
    pass

    ###########################################################

def get_commonbpm(key1, key2, commonbpms):
    i2 = -1
    for i, bpm in enumerate(commonbpms["NAME"]):
        if bpm == key1:
            return i, bpm
        if bpm == key2:
            i2 = i
    if i2 == -1:
        return -1, None
    return i2, commonbpms["NAME"][i2]
