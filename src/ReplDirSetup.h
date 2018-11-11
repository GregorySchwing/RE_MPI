/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ReplDirSetup.h
 * Author: greg
 *
 * Created on November 10, 2018, 7:25 PM
 */

#include <dirent.h>
#include <mpi.h>

#include "repl_ex.h"

#ifndef REPLDIRSETUP_H
#define REPLDIRSETUP_H
class ReplDirSetup{
public:


ReplDirSetup(int temperature, ReplicaExchangeParameters replExParams){
    ReplDirSetup::setupReplicaDirectories(temperature, replExParams);
}

 void setupReplicaDirectories(int temperature, ReplicaExchangeParameters replExParams){   
     
    replica_temp << "temp_" << temperature;
    std::string replica_directory = replica_temp.str();

    mkdirWrapper(replExParams.multisimTitle, replica_directory);
 }
 
 void mkdirWrapper(std::string multisim_directory_name, string replica_directory_name){

    replica_stream << multisim_directory_name << "/" 
                << replica_directory_name;
    std::string replica_directory_path = replica_stream.str();
     
    printf("Creating directory : %s\n", multisim_directory_name.c_str());
    mkdir(multisim_directory_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(replica_directory_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    path_to_replica_directory = replica_directory_path;
    
 }
    
private:
    std::stringstream replica_temp;
    std::stringstream replica_stream;
    std::string path_to_replica_directory;
};
#endif /* REPLDIRSETUP_H */

