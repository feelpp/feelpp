#ifndef _FASTMARCHINGDOFSTATUS_HPP
#define _FASTMARCHINGDOFSTATUS_HPP 1

namespace Feel {

struct FastMarchingDofStatus
{
    typedef uint8_type status_type;
    
    static constexpr status_type FAR = 1 << 0;
    static constexpr status_type CLOSE = 1 << 1;
    static constexpr status_type DONE = 1 << 2;

    static constexpr status_type FIX = 1 << 3;
    static constexpr status_type OLD = 1 << 4;
    static constexpr status_type NEW = 1 << 5;

    static constexpr status_type CLOSE_OLD = CLOSE | OLD;
    static constexpr status_type CLOSE_NEW = CLOSE | NEW;

    static constexpr status_type DONE_FIX = DONE | FIX;
    static constexpr status_type DONE_OLD = DONE | OLD;
    static constexpr status_type DONE_NEW = DONE | NEW;

    FastMarchingDofStatus() = default;
    FastMarchingDofStatus( status_type const& s ): status(s) {}
    operator status_type() const { return status; }

    bool operator&( status_type const s ) const
    {
        return ( this->status & s );
    }
    FastMarchingDofStatus & operator&=( status_type const s )
    {
        this->status &= s;
        return *this;
    }

    bool operator|( status_type const s ) const
    {
        return ( this->status | s );
    }
    FastMarchingDofStatus & operator|=( status_type const s )
    {
        this->status |= s;
        return *this;
    }

    bool operator^( status_type const s ) const
    {
        return ( this->status ^ s );
    }

    FastMarchingDofStatus operator~() const
    {
        return FastMarchingDofStatus( ~( this->status ) );
    }

    status_type status;
};

}

#endif

