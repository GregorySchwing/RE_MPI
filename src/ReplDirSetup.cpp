/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <dirent.h>
#include <mpi.h>

ReplDirSetup::ReplDirSetup(int temperature){
    ReplDirSetup::setupReplicaDirectories(temperature);
}

ReplDirSetup::~ReplDirSetup(){
        delete replica_directory;
        delete directory_stream;
        delete directory_name;
        delete path_stream;
        delete path_string;    
}

 void ReplDirSetup::setupReplicaDirectories(int temperature){   
     
    directory_stream << "temp_" << temperature;
    directory_name = directory_stream.str();
    replica_directory = opendir(directory_name.c_str());
    
    if(replica_directory){
    
        //printf("Directory already exists : %s\n", directory_name.c_str());
        /* Do whatever gromacs does here w the backups #1, #2, ect.
            path_stream << "./" << out.statistics.settings.uniqueStr.val << ".console_out";
            path_string = path_stream.str();
        */
        path_stream << "./" << directory_name << "/replica_log.txt";
        path_string = path_stream.str();
    } else {
        printf("Creating directory : %s\n", directory_name.c_str());
        const int dir_err = mkdir(directory_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err){
            printf("Error creating directory! Writing in pwd\n");
            path_stream << "./" << directory_name << "_replica_log.txt";
            path_string = path_stream.str();
        } else {
            path_stream << "./" << directory_name << "/" << "replica_log.txt";
            path_string = path_stream.str();
        }
    }
 }