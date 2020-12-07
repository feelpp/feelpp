/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/terminalsize.hpp>
#include <feel/feelcore/utility.hpp>

#if defined(FEELPP_OS_WINDOWS)
#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <Windows.h>
#endif
#elif defined(FEELPP_OS_MACOS) || defined(FEELPP_OS_LINUX)
#include <sys/ioctl.h>
#include <fcntl.h>/* open */
#include <unistd.h>/* close */
#endif // Windows/Linux

namespace Feel
{
std::pair<int,int> get_terminal_size()
{
#if defined(FEELPP_OS_WINDOWS)
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    int width = (int)(csbi.dwSize.X);
    int height = (int)(csbi.dwSize.Y);
#elif defined(FEELPP_OS_MACOS)  || defined(FEELPP_OS_LINUX)
    struct winsize w;
#if 0
    //ioctl(fileno(stdout), TIOCGWINSZ, &w);
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w); // not work with mpiexec
#else
    // https://rosettacode.org/wiki/Terminal_control/Dimensions#C
    int fd = open("/dev/tty", O_RDWR);
    if (fd < 0 )
        return std::make_pair(0,0);
    if (ioctl(fd, TIOCGWINSZ, &w) < 0)
        return std::make_pair(0,0);
    close(fd);
#endif
    int width = (int)(w.ws_col);
    int height = (int)(w.ws_row);
#endif // Windows/Linux
    return std::make_pair(width,height);
}

} // Feel
