#!/usr/bin/env bash
# Memory and load snapshot (the WSL VM is shared with another session).
free -m
uptime
ps -eo pid,rss,etime,args --sort=-rss | head -8 | cut -c1-200
