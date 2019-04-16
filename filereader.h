#ifndef FILEREADER_H
#define FILEREADER_H

#include <zlib.h>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <iostream>

namespace util{
    /** Class to hold a file reader in the format .gz or plain text
     */
    class FileReader{
        std::string mFileName;  ///< name of file
        gzFile mGzipFile;       ///< gzFile to store opened mZipped file handler
        FILE* mFile;            ///< FILE pointer to store opened plain file handler
        bool mZipped;           ///< the file is gzipped if true
        char* mBuf;             ///< mBuffer to store a chunk of characters read from mGzipFile or mFile
        int mBufDataLen;        ///< the length of characters read into the mBuffer after the last read
        int mBufUsedLen;        ///< the length of characters already consumed in the mBuffer 
        bool mStdinMode;        ///< read from stdin if true
        bool mNoLineBreakAtEnd; ///< the file has no '\n' as a line break at the last line if true
        size_t mReadBufSize;    ///< the mBuffer size used to read
        
        public:
            /** Construct a file reader with filename
             * @param filename Name of the file
             */
            FileReader(const std::string& filename);
            
            /** FileReader Destructor
             */
            ~FileReader();
           
            /** Tell whether the file is zipped or not
             * @return true if the fasq file is zipped 
             */
            bool isZipped();
    
            /** Get the number of Bytes read from(write to) the file in the meantime, store it in bytesRead
             *  Get the total number of Bytes in the file, store it in bytesTotal
             */
            void getBytes(size_t& bytesRead, size_t& bytesTotal);
    
            /** tell whether the file is a zipped file
             * @param filename file name
             * @return true if the file is a zipped file(with suffix ".gz")
             */
            static bool isZippedFile(const std::string& filename);
            
            /** read one line into line from buffer
             * @param line strin to store result
             * @return true if read successful
             */
            bool getline(std::string& line); 
    
        private:
            /** get just one line from the mBuf, update the mBuf if needed
             */
            std::string getlineFromBuffer();
            
            /** Tell whether the file has no '\n' as a line break at the last line
             * @return true if the file has no '\n' as a line break at the last line
             */ 
            bool hasNoLineBreakAtEnd();

            /** Tell whether the FileReader has reach the endof file
             * @return true if eof reached
             */
            bool eof();
    
            /** initialize the FileReader:
             * 1, open file and store file handler into mGzipFile or mFile or read from stdin
             * 2, set the starting position for the next read on compressed file stream file to the beginning of file 
             * 3, update the file format mZipped
             * 4, call readToBuf() to try to fill the mBuf from first reading
             */ 
            void init();
            
            /** close the FileReader:
             * 1, close file handler
             * 2, set file handler to NULL
             */
            void close();
    
            /** trim \n, \r or \r\n in the tail of the line
             */
            void clearLineBreaks(char* line);
            
            /** try a reading action to fill the mBuf,
             * update mBufDataLen to Bytes read during this reading action
             * reset mBufUsedLen to zero
             * if read the last line(mBuf is not filled), update mNoLineBreakAtEnd
             */
            void readToBuf();
    };
}

#endif
