#ifndef TRANSFORM_FACTORY_HPP
#define TRANSFORM_FACTORY_HPP

#include "../instance/instance.hpp"
#include "../input/input.hpp"

#include "nonstrategic/nonstrategictransform.hpp"
#include "strategic/strategictransform.hpp"
#include "strategic/compact/compacttransform.hpp"
#include "strategic/enum/enumtransform.hpp"
#include "strategic/simple/simpletransform.hpp"

class TransformFactory {
    public:
        inline AbstractTransform* createTransform(const Input & input, const Instance & instance, int pr, Solver *solver){
            Input::Transform chosenTransform = input.getTransform();
            switch (chosenTransform)
            {
                case Input::Transform::NonStrategic:{
                    return new NonStrategicTransform(input,instance,pr,solver);
                    break;
                }
                
                case Input::Transform::Compact :{
                    return new CompactTransform(input,instance,pr,solver);
                    break;         
                }

                case Input::Transform::Enumeration :{
                    return new EnumTransform(input,instance,pr,solver);
                    break;         
                }

                case Input::Transform::Simple :{
                    return new SimpleTransform(input,instance,pr,solver);
                    break;         
                }

                default:{
                    throw std::string("ERROR: Invalid Transform option in command line.");
                    exit(0);
                    break;
                }
            }
            return NULL;
        }
};

#endif