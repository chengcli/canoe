# Linux user guide

This simple manual is intended for a Mac user or a Windows user with WSL (windows
subsustem for linux) installed. Everything is done in a terminal.

## Log on to server
~~~bash
ssh -Y <user_name>@<server_address>
~~~
You will be prompted to enter your password

## Change your password
~~~bash
passwd
~~~
You are provided with an initial password but you should change your password when you
log on for the first time. After you execute the command, you will be asked to enter
your old password first and then the new password.

## Navigate the directories
### Change directory
~~~bash
cd <dir>
~~~

### Change to parent directory
~~~bash
cd ..
~~~

### List files in the directory
~~~bash
ls
~~~
`ls` comes with many options. For example,
~~~bash
ls -l
~~~
will list files vertically and
~~~bash
ls -al
~~~
will list all files including hiden ones

### copy a file
~~~bash
cp <old_file> <new_file>
~~~

### move a file
~~~bash
mv <old_file> <new_file>
~~~

### remove a file
~~~bash
rm <old_file>
~~~

### make a soft link
~~~bash
ln -s <old_file> <new_file>
~~~
A soft link is a pointer to an existing file. It allows your to access the target file without
making a copy. It is better to use absolute path in creating a soft link.

## Edit a file
The basic text editor of a linux system is **vi**. Use the following command to open a
file and edit:
~~~bash
vi <file_name>
~~~
The **vi** editor has two modes: *command mode* and *insert mode*.
### Command mode
The **vi** editor opens in the *command mode* by default, and it only understands commands like
move the cursor, cut, copy and paste the text.
This mode also saves the changes you have made to the file or quits the editor.

### Insert mode
The *insert mode* is for inserting text in the file. You can switch to the *insert mode* from
the *command mode* by pressing `i` on the keyboard. Once you are in the insert mode, any
key would be taken as an input for the file on which you are currently working. To return
to the *command mode*, you need to press the `Esc` key.

### Save changes
When you are in the *command mode*, type
~~~bash
:w
~~~
to save the changes.

### Quit **vi**
When you are in the *command mode*, type
~~~bash
:q
~~~
to quit the **vi** editor.
Combining the two previous commands, you can use
~~~bash
:wq
~~~
to save and quit.

### Quit without saving
To quit without saving, go to the *command mode* and type
~~~bash
:q!
~~~
This will let you leave **vi** without saving the changes.

## Transfer files
On your local computer, use the following command for transfering files on the linux
server to your local computer:
~~~bash
scp <user_name>@<server_address>:<absolute_path_to_file> ./
~~~
In the linux system, `.` represents the current direcotry.
Similarly, you can transfer your local file to a linux server by executing:
~~~bash
scp <local_file> <user_name>@<server_address>:~/
~~~
where `~` represents your home directory on the linux server.
