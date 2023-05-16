import logging
from sys import stdout


def init_logger(verbose: bool = False):
    """
    Initialises a logging object, accessible using the __name__ variable.
    """

    # Decide on logging level
    if verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s.%(msecs)03d %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%m-%d %H:%M:%S",
    )


class ProgressBar:
    """
    Object for displaying and udating a progress bar in the terminal. Example of
    outputs:

    Progress: [..................................................] 0.0 % (0/100)
    Progress: [=========================.........................] 50.0 % (50/100)
    Progress: [==================================================] 100.0 % (100/100)
    """

    def __init__(self, maximum: int, bar_size: int = 50, counter: int = 1):
        """Constructor

        :param maximum: End-value progress is expected to reach
        :type maximum: int
        :param bar_size: Width of displayed bar in terminal, defaults to 50
        :type bar_size: int, optional
        :param counter: Integer value with which to increment, defaults to 1
        :type counter: int, optional
        """
        self.maximum = maximum
        self.bar_size = bar_size  # Reduce for narrow terminal window
        self.counter = counter  # Incremented on .update() call

        # Setup
        init_str = f"Progress: [{'.' * bar_size}] 0.0 % (0/{self.maximum})\r"
        stdout.write(init_str)  # Display bar @ zero progress
        stdout.flush()  # Clear, ready for update
        stdout.write("\b" * (len(init_str)))  # Return to start of line
        # pass

    def update(self) -> None:
        """
        Updates a progress bar object by incrementing progress by ++1. Maximum progress
        for progress bar's size is defined on instantiation of new progress_bar object
        so counter does not have to be parsed.
        """
        # Update information
        prog = int(self.bar_size * self.counter / self.maximum)
        pc_complete = round((self.counter / self.maximum) * 100, 1)
        fraction = f"{self.counter}/{self.maximum}"
        not_prog = "." * (self.bar_size - prog)
        update_str = f"Progress: [{'='*prog}{not_prog}] {pc_complete} % ({fraction})\r"

        # Update progress bar
        stdout.write(update_str)

        # Continue if progress < 100%
        if self.counter < self.maximum:  # bar knows there's more to do
            stdout.flush()
            self.counter += 1
        else:  # action has completed...
            stdout.write("\n")  # ... terminate bar
