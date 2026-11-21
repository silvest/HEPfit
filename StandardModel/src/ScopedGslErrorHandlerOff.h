#ifndef HEPFIT_SCOPED_GSL_ERROR_HANDLER_H
#define HEPFIT_SCOPED_GSL_ERROR_HANDLER_H

#include <gsl/gsl_errno.h>
#include <mutex>

// Keep the process-wide GSL handler off until the last integration leaves
// its scope, including early returns and exceptions.
class ScopedGslErrorHandlerOff {
    struct State {
        std::mutex mutex;
        unsigned int users = 0;
        gsl_error_handler_t* restore = nullptr;
    };
    static State& state() {
        static State shared;
        return shared;
    }
public:
    ScopedGslErrorHandlerOff() {
        State& shared = state();
        std::lock_guard<std::mutex> lock(shared.mutex);
        if (shared.users++ == 0) shared.restore = gsl_set_error_handler_off();
    }
    ~ScopedGslErrorHandlerOff() {
        State& shared = state();
        std::lock_guard<std::mutex> lock(shared.mutex);
        if (--shared.users == 0) gsl_set_error_handler(shared.restore);
    }
    // A persistent handler change takes effect when the last temporary scope exits.
    static void setDefaultHandler(gsl_error_handler_t* handler) {
        State& shared = state();
        std::lock_guard<std::mutex> lock(shared.mutex);
        if (shared.users == 0) gsl_set_error_handler(handler);
        else shared.restore = handler;
    }
    ScopedGslErrorHandlerOff(const ScopedGslErrorHandlerOff&) = delete;
    ScopedGslErrorHandlerOff& operator=(const ScopedGslErrorHandlerOff&) = delete;
};

#endif
