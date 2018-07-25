class AccExcitationMode(object):
    FREE, ACD, ADT = range(3)


class Accelerator(object):
    """
    Abstract class to serve as an interface to implement the
    rest of the accelerators.
    """

    # Class methods ###########################################

    @staticmethod
    def get_class_parameters():
        """
        This method should return the parameter list of arguments needed
        to create the class.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def get_class(cls, *args, **kwargs):
        """
        This method should return the accelerator class defined
        in the arguments.
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
    def get_correctors_variables(cls, frm=None, to=None, classes=None):
        """
        Returns the set of corrector variables between frm and to, with classes
        in classes. None means select all.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def get_element_types_mask(cls, list_of_elements, types):
        """
        Return boolean mask for elements in list_of_elements that belong
        to any of the specified types.
        Needs to handle: "bpm", "magnet", "arc_bpm"

        Args:
            list_of_elements: List of elements
            types: Kinds of elements to look for

        Returns:
            Boolean array of elements of specified kinds.

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
        
    # For GetLLM #############################################################
    
    def get_exciter_bpm(self, plane, distance):
        """
        Returns the BPM next to the exciter.
        The accelerator instance knows already which excitation method is used.
        distance: 1=nearest bpm 2=next to nearest bpm
        """
        Accelerator.__raise_definition_error()
        
    def get_important_phase_advances(self):
        return []
    
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
        raise AttributeError()

    def get_elements_tfs(self):
        """
        Returns the elements tfs file.
        """
        Accelerator.__raise_definition_error()
        
    def get_s_first_BPM(self):
        """
        Returns the position of the first BPM in turn by turn acquisition.
        """
        Accelerator.__raise_definition_error()

    def get_k_first_BPM(self, list_of_bpms):
        """
        Returns the position of something in list_of_bpms TODO: ASK ANDREAS
        """
        Accelerator.__raise_definition_error()

    def get_errordefspath(self):
        Accelerator.__raise_definition_error()

    def set_errordefspath(self, path):
        # TODO: Jaime, are there virtual members for python base classes?
        Accelerator.__raise_definition_error()

    # Templates ##############################################################

    @classmethod
    def get_nominal_tmpl(cls):
        """
        Returns template for nominal model (Model Creator)
        """
        Accelerator.__raise_definition_error()

    # LHC only so far: (put it here because mentioned in creator.py)
    # @classmethod
    # def get_best_knowledge_tmpl(cls):
    #     """
    #     Returns template for best knowledge model
    #     """
    #     Accelerator.__raise_definition_error()
    #
    # @classmethod
    # def get_coupling_tmpl(cls):
    #     """
    #     Returns template for model for coupling correction
    #     """
    #     Accelerator.__raise_definition_error()
    #
    # @classmethod
    # def get_segment_tmpl(cls):
    #     """
    #     Returns template for segment model
    #     """
    #     Accelerator.__raise_definition_error()

    @classmethod
    def get_iteration_tmpl(cls):
        """
        Returns template to create fullresponse.
        TODO: only in _prepare_fullresponse in creator! Needs to be replaced by get_basic_seq
        """
        Accelerator.__raise_definition_error()

    # Jobs ###################################################################

    def get_update_correction_job(self, tiwss_out_path, corrections_file_path):
        """
        Returns job (string) to create an updated model from changeparameters input
        (used in iterative correction).
        """
        Accelerator.__raise_definition_error()

    def get_basic_seq_job(self):
        """
        Returns job (string) to create the basic accelerator sequence.
        """
        Accelerator.__raise_definition_error()

    def get_multi_dpp_job(self, dpp_list):
        """
        Returns job (string) for model with multiple dp/p values (in W-Analysis)
        """
        Accelerator.__raise_definition_error()

    ##########################################################################


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
    """ Returns index and bpm name of key1, if found. Otherwise of key2.
    Raises KeyError if none is found.
    TODO: @Andreas, does that need to be here? Also why not use pythonic ways?

    try:
        return list(commonbpms).index(key1), key1
    except ValueError:
        return list(commonbpms).index(key2), key2

    (try-except around the second return for raising KeyError with message)
    """
    i2 = -1
    for i, bpm in enumerate(commonbpms.index):
        if bpm == key1:
            return i, bpm
        if bpm == key2:
            i2 = i
    if i2 == -1:
        raise KeyError("Giving up looking for the closest exciter BPM")
    return i2, commonbpms.index[i2]

