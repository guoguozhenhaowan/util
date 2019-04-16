#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include "util.h"

namespace util{
    /** Class to write to gz file or ofstream */
    class FileWriter{
            std::string mFilename;  ///< output filename
            gzFile mGzFile;        ///< gzFile file handler
            std::ofstream* mStream; ///< pointer to ofstream
            bool mZipped;           ///< output file is mZipped or not
            int mCompressLevel;     ///< compression level for gz file
            bool mNeedClose;        ///< needed to be closed or not
    
        public:
            /** FileWriter constructor
             * @param filename output filename
             * @param compression compression level for gzFile
             */
            FileWriter(const std::string& filename, const int& compression = 3);
            
            /** FileWriter constructor
             * @param mStream pointer to ofstream
             */
            FileWriter(std::ofstream* mStream);
            
            /** FileWriter constructor
             * @param gzfile gzFile handler
             */
            FileWriter(gzFile gzfile);
            
            /** FileWriter destructor
             */
            ~FileWriter();
    
            /** Whether the output file is mZipped or not
             * @return true if output filename ends with .gz
             */
            bool isZipped();
            
            /** write a string to file without append \\n
             * @param str string to be written 
             * @return true if successfully written
             */
            bool writeString(const std::string& str);
            
            /** write a string to file with additional \\n appended 
             * @param linestr string to be written
             * @return true if successfully written
             */
            bool writeLine(const std::string& linestr);
            
            /** write a C string to file
             * @param cstr pointer to char(C string)
             * @param size length of the C string to be written
             * @return true if successfully written
             */
            bool write(char* cstr, size_t size);
            
            /** get filename of output file
             * @return output filename
             */
            std::string getFilename();
    
            /** initialize FileWriter, detect file format and open file handler
             */
            void init();
            
            /** flush buffer and close file handler
             */
            void close();
    };
}
#endif
