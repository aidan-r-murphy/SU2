/*!
 * \file CNeuralNetwork.cpp
 * \brief Implementation of the NeuralNetwork class to be used
 *      for evaluation of multi-layer perceptrons.
 * \author E. Bunschoten
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */
#include "../../../include/toolboxes/multilayer_perceptron/CNeuralNetwork.hpp"
#include "../../../include/toolboxes/multilayer_perceptron/CLayer.hpp"
#include <iostream>
#include "../../../include/toolboxes/multilayer_perceptron/CReadNeuralNetwork.hpp"
using namespace std;


void MLPToolbox::CNeuralNetwork::predict(vector<su2double> &inputs){
    su2double x, x_norm, y, y_norm;
    size_t iNeuron, iLayer, nNeurons_current;
    bool same_point{true};
    for(iNeuron=0; iNeuron<inputLayer->getNNeurons(); iNeuron++){
        x_norm = (inputs[iNeuron] - input_norm[iNeuron].first)/(input_norm[iNeuron].second - input_norm[iNeuron].first);
        if(abs(x_norm - inputLayer->getOutput(iNeuron)) > 0) same_point = false;
        inputLayer->setOutput(iNeuron, x_norm);
    }
    if(!same_point){
        su2double alpha=1.67326324;
        su2double lambda=1.05070098;
        for(iLayer=1; iLayer<n_hidden_layers + 2; iLayer++){
            nNeurons_current = total_layers[iLayer]->getNNeurons();

            for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                x = ComputeX(iLayer, iNeuron);
                total_layers[iLayer]->setInput(iNeuron, x);
            }
            switch (activation_function_types[iLayer])
                {
                case ENUM_ACTIVATION_FUNCTION::SMOOTH_SLOPE:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);
                        if(x > 0){
                            y = x;
                        }else{
                            y = tanh(x);
                        }
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::ELU:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);
                        if(x > 0){
                            y = x;
                        }else{
                            y = exp(x) - 1;
                        }
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::LINEAR:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);
                        y = x;
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::RELU:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);
                        if(x > 0){
                            y = x;
                        }else{
                            y = 0.0;
                        }
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::SWISH:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);
                        y = x / (1 + exp(-x));
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::TANH:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);
                        y = tanh(x);
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::SIGMOID:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);
                        y = 1.0/(1 + exp(-x));
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::SELU:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);
                        if(x > 0.0){
                            y = lambda * x;
                        }else{
                            y = lambda * alpha * (exp(x) - 1);
                        }
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::GELU:
                    
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        x = total_layers[iLayer]->getInput(iNeuron);

                        y = 0.5 * x * (1 + tanh(sqrt(2 / PI_NUMBER) * (x + 0.044715 * pow(x, 3))));
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                case ENUM_ACTIVATION_FUNCTION::NONE:
                    for(iNeuron=0; iNeuron<nNeurons_current; iNeuron++){
                        y = 0.0;
                        total_layers[iLayer]->setOutput(iNeuron, y);
                    }
                    break;
                default:
                    break;
                }
                
        }
    }
    for(iNeuron=0; iNeuron<outputLayer->getNNeurons(); iNeuron++){
        y_norm = outputLayer->getOutput(iNeuron);
        y = y_norm*(output_norm[iNeuron].second - output_norm[iNeuron].first) + output_norm[iNeuron].first;
        
        // Setting output value
        ANN_outputs[iNeuron] = y;
    }
}

MLPToolbox::CNeuralNetwork::CNeuralNetwork(){
    inputLayer = nullptr;
    outputLayer = nullptr;
    n_hidden_layers = 0;
}


void MLPToolbox::CNeuralNetwork::defineInputLayer(unsigned long n_neurons){
    inputLayer = new CLayer(n_neurons);
    inputLayer->setInput(true);
    input_norm.resize(n_neurons);
}

void MLPToolbox::CNeuralNetwork::defineOutputLayer(unsigned long n_neurons){
    outputLayer = new CLayer(n_neurons);
    output_norm.resize(n_neurons);
}

void MLPToolbox::CNeuralNetwork::push_hidden_layer(unsigned long n_neurons){
    CLayer *newLayer = new CLayer(n_neurons);
    hiddenLayers.push_back(newLayer);
    n_hidden_layers ++;
}

void MLPToolbox::CNeuralNetwork::sizeWeights(){
    
    total_layers.resize(n_hidden_layers + 2);
    total_layers.at(0) = inputLayer;
    for(size_t iLayer=0; iLayer<n_hidden_layers; iLayer++){
        total_layers.at(iLayer + 1) = hiddenLayers.at(iLayer);
    }
    total_layers.at(total_layers.size()-1) = outputLayer;

    weights_mat.resize(n_hidden_layers+1);
    weights_mat[0].resize(hiddenLayers[0]->getNNeurons(), inputLayer->getNNeurons());
    for(size_t iLayer=1; iLayer<n_hidden_layers; iLayer++){
        weights_mat[iLayer].resize(hiddenLayers[iLayer]->getNNeurons(), hiddenLayers[iLayer-1]->getNNeurons());
    }
    weights_mat[n_hidden_layers].resize(outputLayer->getNNeurons(), hiddenLayers[n_hidden_layers-1]->getNNeurons());

    ANN_outputs = new su2double[outputLayer->getNNeurons()];
}   

void MLPToolbox::CNeuralNetwork::setWeight(unsigned long i_layer, unsigned long i_neuron, unsigned long j_neuron, su2double value){
    weights_mat[i_layer][j_neuron][i_neuron] = value;
}
void MLPToolbox::CNeuralNetwork::displayNetwork(){
    /*!
    TODO:
    find way of displaying network architectures
    */
}

void MLPToolbox::CNeuralNetwork::setActivationFunction(unsigned long i_layer, string input)
{
    if(input.compare("linear") == 0){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::LINEAR;
        return;
    }
    if(input.compare("relu") == 0){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::RELU;
        return;
    }
    if((input.compare("smooth_slope") == 0)){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::SMOOTH_SLOPE;
        return;
    }
    if((input.compare("elu") == 0)){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::ELU;
        return;
    }
    if(input.compare("swish") == 0){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::SWISH;
        return;
    }
    if(input.compare("sigmoid") == 0){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::SIGMOID;
        return;
    }
    if(input.compare("tanh") == 0){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::TANH;
        return;
    }
    if(input.compare("selu") == 0){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::SELU;
        return;
    }
    if(input.compare("gelu") == 0){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::GELU;
        return;
    }
    if(input.compare("none") == 0){
        activation_function_types[i_layer] = ENUM_ACTIVATION_FUNCTION::NONE;
        return;
    }
    SU2_MPI::Error("Unknown activation function type: " + input, CURRENT_FUNCTION);
    return;
}