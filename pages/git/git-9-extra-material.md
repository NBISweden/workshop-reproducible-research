The following extra material contains some more advanced things you can do with
git and the command line in general, which is not part of the main course
materials. All the essential skills of git are covered by the previous
sections; the material here should be considered tips and tricks from people
who use git every day. You thus don't need to use these things unless you want
to, and you can even skip this part of the lesson if you like!

## Decorating your prompt

When you are working on the command line interface (CLI), you will usually have
some small pieces of information relating to your current directory, the name
of the computer or host you're working on, and so forth. You've probably
already seen your prompt while working with git throughout this lesson, but
here's an example of what one might look like:

```no-highlight
erikfmbp:~/teaching/workshop-reproducible-research erik.fasterius $
```

The above prompt contains the name of the computer, a colon, the current
working directory, the username and a dollar-sign; it is stored in the
variable `PS1`. You can type `echo $PS1` to see what variables your prompt
is made up of; the above example contains `\h:\W \u\$`, where `\h` is the
hostname, `\W` the working directory and `\u` the username.

Some people like to also show the current branch on their prompt, thus avoiding
having to type `git branch` continuously. There are several ways you might do
this, and we're only presenting one of them here: a bash function.

```bash
git_branch() {
     git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'
}
```

This function does a number of things:

1. Ejects the error message from git if the current directory isn't a part of a
   git repository into `/dev/null` (_i.e._ into nothing).
2. Find the current branch by searching for a line that starts with `*` (*i.e.*
   the current branch) using the command line program `sed`.
3. Put the current branch into parentheses with a space before it.

We can then build our new prompt by adding this function into it:

```bash
# The first part of the old prompt
PS1='\h:\W \u'

# Add the git branch
PS1=$PS1'$(git_branch)'

# Add the last part of the old prompt
PS1=$PS1' \$'
```

Now you should see the current git branch on your prompt! The only problem now
is that this only works for your current session: once you restart your CLI
you'll have to re-define your prompt again. This can be circumvented, though.
What you need to do is to add the code defining your prompt into your so-called
bash profile: `~/.bash_profile`. Every time you load a new CLI session this
file is read and any code inside it is executed. You might already have this
file, so make sure you don't overwrite it!

## Bash aliases for git

Some git commands are used over and over again when working with git, such as
`git status`. Some people like to have aliases (*i.e.* shortcuts) for these
common commands. Here is a small list of such aliases that you may find useful
or, even better, might inspire you to create your own! Add them to your
`~/.bash_profile` as above, so that they're available across sessions.

```bash
# Basic git commands
alias gb='git branch'
alias ga='git add'
alias gd='git diff'
alias gcm='git commit'
alias gp='git push'
alias gu='git pull'
alias gm='git merge'
alias gco='git checkout'
alias gl='git log'

# Git status in short format
alias gst='git status -s'

# Add and commit all tracked and modified files
alias gca='git commit -a'

# Create and checkout a new branch
alias gcob='git checkout -b'

# Git log with one line per commit
alias glo='git log --oneline'
```

## Remote connections with SSH

Throughout this tutorial we have been using normal HTTPS to connect with remote
repositories, which means you need to provide your username and password. If you
work with remote repositories a lot you might want to be able to skip writing in
your user credentials all the time, which is something you can do with SSH
(Secure Shell) keys. This basically entails setting up a pair of keys: one
private and one public. You keep the private key on your local computer and give
the public key to anywhere you want to be able to connect to, *e.g.* GitHub. The
public key can be used to encrypt messages that *only* the corresponding private
key can decrypt. A simplified description of how SSH authentication works goes
like this:

1. The client (*i.e.* the local computer) sends the ID of the SSH key pair it
   would like to use for authentication to the server (*e.g.* GitHub)
2. If that ID is found, the server generates a random number and encrypts this
   with the public key and sends it back to the client
3. The client decrypts the random number with the private key and sends it back
   to the server

Notice that the private key always remains on the client's side and is never
transferred over the connection; the ability to decrypt messages encrypted with
the public key is enough to ascertain the client's authenticity. This is in
contrast with using passwords, which are themselves send across a connection
(albeit encrypted). It is also important to note that even though the keys come
in pairs it is impossible to derive the private key from the public key.

So, how do you setup SSH authentication between your local computer and GitHub?
GitHub actually has excellent and detailed instructions for exactly this at the
[GitHub website](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh),
so we won't repeat that here. Simply follow those instructions and you'll be
good to go! If you want to read more details about how SSH authentication work
you can check out [this website](https://www.digitalocean.com/community/tutorials/understanding-the-ssh-encryption-and-connection-process),
which has more in-depth information than we provide here.

## Forking

When you want to work on an Open Source project that is available on *e.g.*
GitHub, you usually don't have permission to directly push code to the project's
repository - this is so that the project's maintainers are the only ones that
can directly change anything in their codebase. How do you then contribute to
projects that don't allow you to push your code to their repository? Simple: use
*forking*!

Forking is when you make your own copy of a repository on your GitHub account,
which you will then have permissions to change as you see fit. You can then
create pull requests from your fork to the original repository, rather than
pushing code to a new branch and making a pull request from that. Working with
forks just adds an additional step to the whole workflow: instead of being
"clone; code and commit changes on a new branch; push branch to remote; pull
request from branch" it becomes "fork; clone; code and commit changes; push code
to fork; pull request from fork".

You might also want to do a fork of a project simply because you want to have
your own copy of it as well, without ever having the intention of changing it.
This is, of course, perfectly fine as well, but do keep in mind that developers
are usually quite happy to incorporate new changes from contributors if they are
reasonable and fulfil a purpose and add functionality to the project. It is
quite common that you have a use-case the maintainer didn't think of before, and
that you've helped the project grow by contributing your code!
