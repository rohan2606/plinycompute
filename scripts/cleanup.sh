#!/usr/bin/env bash
#  Copyright 2018 Rice University
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#  ========================================================================    
set -o errexit 

usage() {
    echo ""    
    echo -e "\033[33;31m""    "Warning: This script deletes stored data, use with care!"\e[0m"

    cat <<EOM

    Description: This script deletes all PlinyCompute storage, catalog metadata,
    and kills both pdb-manager and pdb-worker processes in one machine if 
    cluster-type=standalone or the entire cluster if cluster_type=distributed.

    Usage: scripts/$(basename $0) param1 param2

           param1: <pem_file>
                      Specify the private key to connect to other machines in 
                      the cluster; the default is conf/pdb-key.pem
           param2: <cluster_type>
                      Specify the type of cluster; {'standalone', 'distributed'}

EOM
   exit -1;
}

[ -z $1 ] && [ -z $2 ] && { usage; } || [[ "$@" = *--help ]] && { usage; } || [[ "$@" = *-h ]] && { usage; }

pem_file=$1
cluster_type=$2
user=ubuntu
ip_len_valid=3
pdb_dir=$PDB_INSTALL
testSSHTimeout=3

if [ ! -f ${pem_file} ]; then
    echo -e "Pem file ""\033[33;31m""'$pem_file'""\e[0m"" not found, make sure the path and file name are correct!"
    exit -1;
fi

scripts/cleanupNode.sh

# By default disable strict host key checking
if [ "$PDB_SSH_OPTS" = "" ]; then
  PDB_SSH_OPTS="-o StrictHostKeyChecking=no"
fi

if [ -z ${pem_file} ];
then
  PDB_SSH_OPTS=$PDB_SSH_OPTS
else
  PDB_SSH_OPTS="-i ${pem_file} $PDB_SSH_OPTS"
fi

if [ "$cluster_type" != "standalone" ] && [ "$cluster_type" != "distributed" ];
   then echo "ERROR: the value of cluster_type can only be either: 'standalone' or 'distributed'";
   exit -1;
fi

# parses conf/serverlist file
if [ "$cluster_type" = "standalone" ];then
   conf_file="conf/serverlist.test"
else
   conf_file="conf/serverlist"
fi

while read line
do
   [[ $line == *#* ]] && continue # skips commented lines
   [[ ! -z "${line// }" ]] && arr[i++]=$line # include only non-empty lines
done < $PDB_HOME/$conf_file

length=${#arr[@]}
echo "There are $length servers defined in $PDB_HOME/$conf_file"

for (( i=0 ; i<=$length ; i++ ))
do
   ip_addr=${arr[i]}
   if [ ${#ip_addr} -gt "$ip_len_valid" ];then
      only_ip=${ip_addr%:*}
      # checks that ssh to a node is possible, times out after 3 seconds
      if [[ ${ip_addr} != *":"* ]];then
         nc -zw$testSSHTimeout ${ip_addr} 22
      else
         nc -zw$testSSHTimeout ${only_ip} 22
      fi
      if [ $? -eq 0 ];then
            echo -e "\n+++++++++++ cleanup server: $ip_addr"
         if [[ ${ip_addr} != *":"* ]];then
            ssh $PDB_SSH_OPTS $user@$ip_addr "cd $pdb_dir; scripts/cleanupNode.sh"
         else
            ./scripts/cleanupNode.sh
         fi
      else
         echo "Cannot clean server with IP address: ${ip_addr}, connection timed out on port 22 after $testSSHTimeout seconds."
      fi
    fi
done


