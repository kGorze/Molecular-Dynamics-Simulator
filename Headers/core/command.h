//WRITTEN BY: github.com/kgorze

#ifndef COMMAND_H
#define COMMAND_H



class command {
public:
    virtual ~command() {}
    virtual void execute() = 0;
};



#endif //COMMAND_H
