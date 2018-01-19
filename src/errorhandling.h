/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Error handling macros.
 *
 * @author Carlos Fernandes
 * @author Nuno Fachada
 */

/**
 * Macro for terminating program with error condition.
 *
 * @param[in] err_msg Error message.
 * @param[in] ... Message parameters.
 */
#define ERROR_EXIT(err_msg, ...) \
	do { \
		fprintf(stderr, err_msg "\n", __VA_ARGS__); \
		exit(EXIT_FAILURE); \
	} while (0)
