#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>
#include <zip.h>
#include <fmt/core.h>
#include <feel/feelcore/feel.hpp>
namespace Feel
{
    bool extractZipFile( const std::string& zipFilePath, const std::string& extractionDir );
    void cleanupTemporaryDirectory( const std::string& extractionDir );
}