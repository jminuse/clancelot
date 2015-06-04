import os, sys

import smtplib
from email.mime.text import MIMEText

msg = MIMEText(sys.argv[1] + " has finished")

msg['Subject'] = "Simulation Notification"
msg['From'] = "hherbol@gmail.com"
msg['To'] = sys.argv[2]

# Send the message via our own SMTP server, but don't include the
# envelope header.
s = smtplib.SMTP('localhost')
s.sendmail("hherbol@gmail.com", [sys.argv[2]], msg.as_string())
s.quit()