# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date:   2017-12-05 12:01:21
# @Last modified by:   Brian Cherinka
# @Last Modified time: 2017-12-05 12:19:32

from __future__ import print_function, division, absolute_import


class PlateschedulerError(Exception):
    """A custom core Platescheduler exception"""

    def __init__(self, message=None):

        message = 'There has been an error' \
            if not message else message

        super(PlateschedulerError, self).__init__(message)


class PlateschedulerNotImplemented(PlateschedulerError):
    """A custom exception for not yet implemented features."""

    def __init__(self, message=None):

        message = 'This feature is not implemented yet.' \
            if not message else message

        super(PlateschedulerNotImplemented, self).__init__(message)


class PlateschedulerAPIError(PlateschedulerError):
    """A custom exception for API errors"""

    def __init__(self, message=None):
        if not message:
            message = 'Error with Http Response from Platescheduler API'
        else:
            message = 'Http response error from Platescheduler API. {0}'.format(message)

        super(PlateschedulerAPIError, self).__init__(message)


class PlateschedulerApiAuthError(PlateschedulerAPIError):
    """A custom exception for API authentication errors"""
    pass


class PlateschedulerMissingDependency(PlateschedulerError):
    """A custom exception for missing dependencies."""
    pass


class PlateschedulerWarning(Warning):
    """Base warning for Platescheduler."""


class PlateschedulerUserWarning(UserWarning, PlateschedulerWarning):
    """The primary warning class."""
    pass


class PlateschedulerSkippedTestWarning(PlateschedulerUserWarning):
    """A warning for when a test is skipped."""
    pass


class PlateschedulerDeprecationWarning(PlateschedulerUserWarning):
    """A warning for deprecated features."""
    pass
