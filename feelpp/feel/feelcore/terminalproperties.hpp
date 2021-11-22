/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#pragma once

#include <optional>

namespace Feel
{

/**
 * @brief describe terminal properties if available : width, height
 */
class TerminalProperties
{
public :
    TerminalProperties();
    static bool hasWidth() { return S_width && *S_width > 0; }
    static bool hasHeight() { return S_height && *S_height > 0; }
    static int width() { return *S_width; }
    static int height() { return *S_height; }

    static void updateTerminalSize( bool updateEvenIfAlreadyDefined = false );
private :
    inline static std::optional<int> S_width;
    inline static std::optional<int> S_height;
};

} // Feel
