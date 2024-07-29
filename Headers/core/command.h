//
// Created by konrad_guest on 18/05/2024.
//

#ifndef COMMAND_H
#define COMMAND_H



class command {
public:
    virtual ~command() {}
    virtual void execute() = 0;
};



#endif //COMMAND_H
