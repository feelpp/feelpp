#pragma once

namespace Feel
{
    bool extractZipFile( const std::string& zipFilePath, const std::string& extractionDir );
    void cleanupTemporaryDirectory( const std::string& extractionDir );
}