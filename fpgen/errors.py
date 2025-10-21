"""Project error types exposed to callers."""

class WignerFPError(ValueError):
    """Raised for domain/representation errors in the Wigner FP toolkit."""
    pass

class ParseError(ValueError):
    """Raised for expression parsing/tokenization errors."""
    pass