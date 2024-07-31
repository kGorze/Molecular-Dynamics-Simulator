//BASED ON THE GROMACS SOURCE CODE github.com/gromacs/gromacs


#ifndef MDS_UTILITY_ERRORCODES_H
#define MDS_UTILITY_ERRORCODES_H

namespace mds{
enum ErrorCode
{
eeOK,
eeOutOfMemory,
eeFileNotFound,
eeFileIO,
eeInvalidInput,
eeInconsistentInput,
eeTolerance,
eeInstability,


eeNotImplemented,
eeInvalidValue,
eeInvalidCall,
eeInternalError,
eeAPIError,
eeRange,
eeParallelConsistency,
eeModularSimulator,
eeUnknownError,
};

const char* getErrorCodeString(int errorcode);


}
#endif //MDS_UTILITY_ERRORCODES_H
