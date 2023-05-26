/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <exception>
#include <string>

class AnnotationException : public std::logic_error
{
public:
    AnnotationException(const std::string& s) : std::logic_error(s) {}
};