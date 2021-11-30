#!/bin/bash

set -e

source ~/jenkins-env.bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 4 ]; then
  echo "Start pod for tool installations"
  echo "Usage: $0 JOB_NAME BUILD_NUNMBER IMAGE BUNDLE_COLLECTION_VERSION"
  exit 0
fi

JOB_NAME="$1"
BUILD_NUMBER="$2"
image="$3"
BUNDLE_COLLECTION_VERSION="$4"

adjusted_job_name="$(echo "$JOB_NAME" | tr '_' '-' )"
name="tool-install-$adjusted_job_name-$BUILD_NUMBER"

# umount old tools-bin if exists
if grep -qs "/mnt/data/$name/tools-bin " /proc/mouns; then
  echo "umount tools-bin"
  sudo umount /mnt/data/$name/tools-bin
fi
    
# create an upper dir for the union file system"
sudo mkdir -p /mnt/data/$name/tools-bin-upper
sudo chown ubuntu:ubuntu /mnt/data/$name/tools-bin-upper

echo "mount tools-bin: $BUNDLE_COLLECTION_VERSION"
sudo mkdir -p /mnt/data/$name/tools-bin

# using unionfs-fuse, because overlayfs refused to make modifications to subdirectories when used on top of NFS 
unionfs_cmd="unionfs-fuse -o cow -o allow_other /mnt/data/$name/tools-bin-upper=RW"

for f in /mnt/artefacts/collect_bundles/$BUNDLE_COLLECTION_VERSION/*; do
  unionfs_cmd="$unionfs_cmd:$f=RO"
done

unionfs_cmd="$unionfs_cmd /mnt/data/$name/tools-bin"
#echo $unionfs_cmd
sudo bash -c "$unionfs_cmd"

# delete old deployment if exists
if kubectl get deployment $name > /dev/null 2>&1; then
  kubectl delete deployment $name
fi

cat <<EOF | kubectl apply -f -
apiVersion: apps/v1
kind: Deployment
metadata:
  name: $name
spec:
  replicas: 1
  selector:
    matchLabels:
      app: $name
  template:
    metadata:
      labels:
        app: $name
    spec:
      containers:
      - name: $name
        image: docker-registry.rahti.csc.fi/chipster-images-beta/$image
        volumeMounts:
        - name: tools-bin
          mountPath: /mnt/tools
        - name: artefacts
          mountPath: /mnt/artefacts
      volumes:
      - name: tools-bin
        hostPath:
          path: /mnt/data/$name/tools-bin
          type: Directory
      - name: artefacts
        hostPath:
          path: /mnt/artefacts
          type: Directory
EOF

echo "waiting pod to start"
for i in $(seq 30); do
  if kubectl get pod 2> /dev/null | grep $name | grep Running > /dev/null; then     
    echo ""
    break
  fi
  printf "."
  sleep 1
done
