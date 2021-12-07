//
//  MLException.h
//  VAMaximumLikelihood
//
//  Created by Hugh Dickinson on 10/28/14.
//  Copyright (c) 2014 Josh Cardenzana. All rights reserved.
//

#ifndef VAMaximumLikelihood_MLException_h
#define VAMaximumLikelihood_MLException_h

#include <stdexcept>
#include <string>

//typedef std::runtime_error MLException;

class MLException : public std::exception {
public:
    explicit MLException() :
        error_msg_("")
        {}
    
    explicit MLException(const std::string& errMsg) :
        error_msg_(errMsg)
        {}
    
    virtual ~MLException() throw() {}
    
    virtual const char* what() const throw ()
        {return error_msg_.c_str();}
    
    
protected:
    // Error message string
    std::string error_msg_ ;
    
} ;

#endif
