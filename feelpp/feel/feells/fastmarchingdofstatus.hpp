#ifndef _FASTMARCHINGDOFSTATUS_HPP
#define _FASTMARCHINGDOFSTATUS_HPP 1

namespace Feel {

struct FastMarchingDofStatus
{
    typedef uint8_type status_type;
    
    static constexpr status_type FAR = 0;

    static constexpr status_type CLOSE_OLD = 1;
    static constexpr status_type CLOSE_NEW = 2;

    static constexpr status_type DONE_FIX = 3;
    static constexpr status_type DONE_OLD = 4;
    static constexpr status_type DONE_NEW = 5;

    FastMarchingDofStatus() = default;
    FastMarchingDofStatus( status_type const& s ): status(s) {}
    operator status_type() const { return status; }

    status_type status;
};

}

#endif

