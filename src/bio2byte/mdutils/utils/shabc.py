from abc import ABCMeta
import tempfile


class ShellInterface(metaclass=ABCMeta):
    
    def __init__(self, **kwargs):
        self._tmpfiles = []
        
    def tmpfile(self, filename, **tempfile_kwargs):
        """ Uses a filename or creates a NamedTempfile """
        if filename is None:
            tmp = tempfile.NamedTemporaryFile(**tempfile_kwargs)
            self._tmpfiles.append(tmp)
            return tmp.name
        else:
            return filename
