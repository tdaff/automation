"""

Abstract classes that can be plugged in to register job managers.

"""

from abc import abstractmethod
from abc import ABCMeta


class JobManager(object):
    """
    Abstract base class for job manager.
    Implemented managers must be subclasses of JobManager and
    implement all the methods with the documented return values.

    """

    __metaclass__ = ABCMeta

    managers = []

    def __init__(self, ):
        self.register_managers()

    @classmethod
    def register_managers(cls):
        for subclass in cls.__subclasses__():
            cls.managers.append(subclass.name)

    @staticmethod
    @abstractmethod
    def submit():
        """
        Submit a job to be run.

        Valid retun values are:
            jobid<str> - for a queued job;
            True<bool> - job has run and completed;
            False<bool> - job has failed.
        """
        return None

    @staticmethod
    @abstractmethod
    def postrun(jobid=None, jobids=None):
        """
        Submit the script to run itself after the job is complete.

        Valid retun values are:
            jobid<str> - for a queued job;
            True<bool> - job has run and completed;
            False<bool> - job has failed.
        """
        return None

    @staticmethod
    @abstractmethod
    def jobcheck(jobid):
        """
        Query the system to see if the job is still running.

        Valid retun values are:
            True<bool> - job is queued or still running;
            False<bool> - job is no longer in the queue.
        """
        return None

    @staticmethod
    @abstractmethod
    def env(job_type):
        """
        Prepare the running environment with any special conditions.
        
        Method has no return value.
        """
        return None
