# 📟 Keeping Terminal Sessions Alive with `screen`

When you run commands on a remote server (via SSH), they usually stop if you close your laptop, lose internet, or log out. `screen` fixes this by letting you create a detachable session that keeps running in the background. You can disconnect and reconnect later without losing progress.

## 1. ✅ Install `screen`

Most Linux servers already have it. If not:

```bash
sudo apt-get install screen
```

_(You can skip this if it's already installed.)_

## 2. 🛠️ Start a New Screen Session

Run:

```bash
screen -S myproject
```

Here, `myproject` is the name of your screen session.

You'll see a welcome message and a blank terminal — you're now inside a screen session.

You can start your long-running job, for example:

```bash
spades.py -s nanopore.fastq.gz
```

## 3. 🔌 Detach from the Session (Keep It Running)

To disconnect but keep the job running, press:

```
Ctrl + A, then D
```

Think of it as:

- `Ctrl + A` → tell `screen` you're issuing a command
- `D` → detach

Now you're back in your normal terminal, but the `screen` session is still alive in the background.

## 4. 📋 List Existing Screen Sessions

Use:

```bash
screen -ls
```

You will see something like:

```
There is a screen on:
    myproject.pts-0.server-name  (Detached)
1 Socket in /run/screen/S-username.
```

## 5. 🔁 Reattach to a Session

To go back into your running session:

```bash
screen -r myproject
```

Or simply:

```bash
screen -r
```

_(If there's only one session)_

## 6. ❌ Kill a Session (When You're Done)

Inside the session, type:

```bash
exit
```

Or press:

```
Ctrl + D
```

This closes the screen and frees up system resources.

## 🧠 Quick Reference

| Action        | Command / Shortcut          |
|---------------|------------------------------|
| Start         | `screen -S name`             |
| Detach        | `Ctrl + A` then `D`          |
| List sessions | `screen -ls`                 |
| Reattach      | `screen -r name`             |
| Kill/Exit     | `exit` (inside the screen)   |
