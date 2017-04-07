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
        by in the arguments.
        """
        Accelerator.__raise_non_implemented()

    @classmethod
    def get_class_from_args(cls, args):
        """
        This method should return the accelerator class defined
        by in the passed command line like arguments.
        """
        Accelerator.__raise_non_implemented()

    # Instance methods ########################################

    def init_from_args(cls, args):
        """
        Instances an accelerator filling the attributes using the
        command line like arguments passed in args.
        """
        Accelerator.__raise_non_implemented()

    def verify_object(self):
        """
        Verifies that this instance of an accelerator is properly
        instantiated.
        """
        Accelerator.__raise_non_implemented()

    @staticmethod
    def __raise_non_implemented():
        raise NotImplementedError(
            """A function that should have been overwritten
            has been called, check stack trace."""
        )

    ###########################################################
