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
    def get_bpms_around_acd(available_bpm_list, model_twiss, exciter="acd"):
        """
        It will return the two best BMPs to use for exciter specified
        in compensation.
        """
        Accelerator.__raise_definition_error()

    @classmethod
    def get_arc_bpms(list_of_elements):
        """
        It will return a list with the elements of @list_of_elements that
        are BPMs of the arcs.
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

    ###########################################################


class AcceleratorDefinitionError(Exception):
    """
    Raised when an accelerator instance is wrongly used, for
    example by calling a method that should have been overwritten.
    """
    pass
