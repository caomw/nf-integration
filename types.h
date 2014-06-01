#ifndef TYPES_H_
#define TYPES_H_

/*! \brief
 * a list of data types commonly used in R4R
 *
 */
enum class ETYPE {  NA = 0,
                    B1U = 1,
                    C1U = 2, C1S = 3, C1U3 = 102,
                    S2U = 4, S2S = 5,
                    I4S = 6, I4U = 7,
                    F4S = 8, F4S3 = 108,
                    L8S = 9, L8U = 10,
                    D8S = 11, D8S3 = 111,
                    STRING = 19 };

template<typename T> ETYPE GetEType();

#endif /* TYPES_H_ */
