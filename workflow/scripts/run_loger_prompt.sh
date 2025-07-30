#!/usr/bin/env expect

# Set the timeout for Expect to wait for a prompt. Increase if needed.
set timeout 25

# Launch Orthologer
spawn sbin/setup.sh go

expect "first prompt from Orthologer:"
send "\r"
expect "2nd prompt from Orthologer:"
send "\r"
expect "3nd prompt from Orthologer:"
send "path2proteome\r"
expect "4nd prompt from Orthologer:"
send "\r"
expect "5nd prompt from Orthologer:"
send "\r"
expect "6nd prompt from Orthologer:"
send "\r"
expect "7nd prompt from Orthologer:"
send "\r"
expect "8nd prompt from Orthologer:"
send "\r"
expect "9nd prompt from Orthologer:"
send "\r"
expect "10nd prompt from Orthologer:"
send "100\r"
expect "11nd prompt from Orthologer:"
send "500\r"
expect "12nd prompt from Orthologer:"
send "\r"
expect "13nd prompt from Orthologer:"
send "\r"
expect "14nd prompt from Orthologer:"
send "/work/FAC/FBM/DEE/mrobinso/moult/giulia/orthologer_v3/orthologer_3.0.2/BRHCLUS-5.1.6/bin\r"


# Wait for the program to finish
expect eof

