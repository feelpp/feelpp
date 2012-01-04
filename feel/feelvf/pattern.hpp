
#ifndef __PATTERN_H
#define __PATTERN_H 1


namespace Feel
{

    //namespace vf
    //{
enum  Pattern
    {
        DEFAULT   = 1 << 0,
        EXTENDED  = 1 << 1,
        COUPLED   = 1 << 2,
        SYMMETRIC = 1 << 3,
        ZERO      = 1 << 4,
        HAS_NO_BLOCK_PATTERN = 1 << 5
    };

    //} // namespace vf

} // namespace Feel


#endif // __PATTERN_H
