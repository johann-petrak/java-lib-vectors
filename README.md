# java-lib-vectors

A very simple and small-footprint library to collect some vector utilities which 
are useful to handle array of double dense vectors in Java. This has mainly be
written to make it a bit easier to deal with word embeddings represented as 
arrays of double. 

This library implements all functions as instance methods, so an instance of VectorUtils
has to be created first. Depending on how the instance is created, the methods will
differ in certain configuration dependent details, e.g. whether an error should be
thrown for certain sitations where no reasonable result can be calculated or just
NaN should get returned instead.

LICENSE: LGPL 2.1
