# Cleaning up

As mentioned before, Docker tends to consume a lot of disk space. In general,
`docker image rm` is used for removing images and `docker container rm` for
removing containers. Here are some convenient commands for cleaning up.

```bash
# Remove unused images
docker image prune

# Remove stopped containers
docker container prune

# Remove unused volumes (not used here, but included for reference)
docker volume prune

# Stop and remove ALL containers
docker container rm $(docker container ls -a -q)

# Remove ALL images
docker image rm $(docker image ls -a -q)
```
