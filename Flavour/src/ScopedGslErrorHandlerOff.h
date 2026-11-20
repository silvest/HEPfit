#ifndef HEPFIT_SCOPED_GSL_ERROR_HANDLER_H
#define HEPFIT_SCOPED_GSL_ERROR_HANDLER_H

#include <gsl/gsl_errno.h>
#include <mutex>

// GSL has one process-wide handler. Restore it when the last overlapping
// integration leaves its scope, including early returns and exceptions.
class ScopedGslErrorHandlerOff {
    struct State {
        std::mutex mutex;
        unsigned int users = 0;
        gsl_error_handler_t* previous = nullptr;
    };
    static State& state() {
        static State shared;
        return shared;
    }
public:
    ScopedGslErrorHandlerOff() {
        State& shared = state();
        std::lock_guard<std::mutex> lock(shared.mutex);
        if (shared.users++ == 0) shared.previous = gsl_set_error_handler_off();
    }
    ~ScopedGslErrorHandlerOff() {
        State& shared = state();
        std::lock_guard<std::mutex> lock(shared.mutex);
        if (--shared.users == 0) gsl_set_error_handler(shared.previous);
    }
    ScopedGslErrorHandlerOff(const ScopedGslErrorHandlerOff&) = delete;
    ScopedGslErrorHandlerOff& operator=(const ScopedGslErrorHandlerOff&) = delete;
};

#endif
