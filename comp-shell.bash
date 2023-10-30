#!/bin/bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -gt 4 ]; then
  echo "Open shell in container image"
  echo "Usage: $0 [--local-image] [--always-pull] [--no-tools] [IMAGE]"
  exit 0
fi

image_repo="docker-registry.rahti.csc.fi/chipster-images-release/"
image_pull_policy="IfNotPresent"

# --local-image
for arg in "$@"; do
  if [[ "$arg" == "--local-image" ]]; then
    image_pull_policy="Never"
  fi
done

# --always-pull
for arg in "$@"; do
  if [[ "$arg" == "--always-pull" ]]; then
    image_pull_policy="Always"
  fi
done

# --no-tools
for arg in "$@"; do
  if [[ "$arg" == "--no-tools" ]]; then
    no_tools=true
  fi
done

# image
image="comp-20-04-r-deps"

# Use the last argument as the image if it doesn't start with '--'
if [ $# -gt 0 ]; then
  # Get the last argument
  last_argument="${!#}"

  # Check if the last argument doesn't start with '--'
  if [[ ! "$last_argument" == --* ]]; then
    image="$last_argument"
  fi
fi

name="comp-shell"

if ! kubectl get pod | grep $name | grep Running >/dev/null 2>&1; then
  deployment=$(
    cat <<EOF
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
        volumeMounts: []
      volumes: []
EOF
  )

  echo "image: ${image_repo}${image}"
  echo "imagePullPolicy: ${image_pull_policy}"

  patch=".spec.template.spec.containers[0].image = \"${image_repo}${image}\" |
    .spec.template.spec.containers[0].imagePullPolicy = \"${image_pull_policy}\""

  TOOLS_BIN_NAME=$(cat ~/values.yaml | yq e - -o json | jq .toolsBin.version -r)
  TOOLS_BIN_HOST_MOUNT_PATH=$(cat ~/values.yaml | yq e - -o json | jq .toolsBin.hostPath -r)
  TOOLS_BIN_PATH="/opt/chipster/tools"

  if [[ $TOOLS_BIN_NAME != "null" && $no_tools != "true" ]]; then
    patch="$patch |
      .spec.template.spec.containers[0].volumeMounts += [{\"name\": \"tools-bin\", \"readOnly\": false, \"mountPath\": \"$TOOLS_BIN_PATH\"}]"

    if [[ $TOOLS_BIN_HOST_MOUNT_PATH != "null" ]]; then
      echo "mount tools-bin from hostPath $TOOLS_BIN_HOST_MOUNT_PATH/$TOOLS_BIN_NAME to $TOOLS_BIN_PATH"
      patch="$patch |
        .spec.template.spec.volumes += [{\"name\": \"tools-bin\", \"hostPath\": { \"path\": \"$TOOLS_BIN_HOST_MOUNT_PATH/$TOOLS_BIN_NAME\", \"type\": \"Directory\"}}]"
    else
      echo "mount tools-bin from PVC $TOOLS_BIN_NAME to $TOOLS_BIN_PATH"
      patch="$patch |
        .spec.template.spec.volumes += [{\"name\": \"tools-bin\", \"persistentVolumeClaim\": { \"claimName\": \"tools-bin-$TOOLS_BIN_NAME\"}}]"
    fi
  else
    echo "do not mount tools-bin"
  fi

  json=$(echo "$deployment" | yq e - -o=json | jq "$patch")

  echo "$json" | kubectl apply -f -

  echo "waiting pod to start"
  for i in $(seq 30); do
    if kubectl get pod 2>/dev/null | grep $name | grep Running >/dev/null; then
      echo ""
      break
    fi
    printf "."
    sleep 1
  done
fi

echo "comp shell, press Ctrl+D to exit"
kubectl exec -it $(kubectl get pod | grep $name | grep Running | cut -d " " -f 1) -- bash

kubectl delete deployment $name
