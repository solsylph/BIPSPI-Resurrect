from Config import Configuration

class MyException(Exception):
  def __init__(self, msg):
    Exception.__init__(self, msg)
    self.msg= msg
  
class NoAvailableForDownloadPDB(MyException):
  def __init__(self, msg):
    MyException.__init__(self, msg)
  
class NoValidPDBFile(MyException):
  def __init__(self, msg):
    MyException.__init__(self, msg)
    
class BadNumberOfResidues(Configuration, MyException):
  def __init__(self, nResidues, partnerId):
    Configuration.__init__(self)
    self.nResidues = nResidues
    self.partnerId = partnerId
    MyException.__init__(self, "Bad number of residues for partner %s: %d. Number of residues must be %d < nResidues < %d"%(
                                 partnerId, nResidues, self.minNumResiduesPartner , self.maxNumResiduesPartner))

  def __reduce__(self):
    # Pickle round-trips an exception via cls(*args); this exception's __init__
    # takes (nResidues, partnerId), not the single message string MyException
    # stores in self.args. Without this, unpickling in the parent process (when
    # a multiprocessing worker raises it) fails with a TypeError and kills the
    # pool's result-handler thread.
    return (self.__class__, (self.nResidues, self.partnerId))
    
class BadSequence(Configuration, MyException):
  def __init__(self, msg):
    Configuration.__init__(self)
    MyException.__init__(self, msg)
    
