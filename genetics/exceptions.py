class MendelianInferenceException(Exception):
    pass

class ParserException(MendelianInferenceException):
    pass

class InvalidObservation(ParserException):
    pass

class InvalidChild(ParserException):
    pass

class AnalysisException(MendelianInferenceException):
    pass

class InvalidState(AnalysisException):
    pass

class NonMendelianPattern(AnalysisException):
    pass