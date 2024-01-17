'''
This file sends emails if there are failing tests or if a previously failing test passes now.
'''
import os, sys
from store import get_version
import logpage

# Import smtplib for the actual sending function
import smtplib

# Import the email modules we'll need
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

logpage.REPO = sys.argv[1]

# Create the body of the message (a plain-text and an HTML version).
text = "HTML only email, please see https://einsteintoolkit.github.io/carpetx-tests for output"

curr_ver = get_version()
summary=f"./records/version_{curr_ver}/build__2_1_{curr_ver}.log"
baseurl = "https://github.com/einsteintoolkit/carpetx-tests"

data = logpage.create_summary(summary)
status = "All Tests Passed"
if data["Number failed"]!=0:
    status="Some Tests Failed"
# Add a table row for each data field
contents = "\n".join([f"<tr><th>{key}</th><td>{data[key]}</td><tr>" for key in data.keys()])

html = f'''<!doctype html>
    <html lang="en">
		<head></head>
        <body>
            <h3 style="text-align:center"><a href="{baseurl}/tree/gh-pages/records/version_{curr_ver}">Build #{curr_ver}</a></h3>
            <table class="table table-bordered " >
            <caption>Summary</caption>
            {contents}
            </table>
            <table>
            <caption>Commits in Last Push</caption>
            {logpage.gen_commits()}
            </table>
            {logpage.gen_diffs(summary)}
		</body>
	</html>
'''

msg = MIMEMultipart('alternative')
msg['Subject'] = f"CarpetX test report: {status}"
msg['From'] = "jenkins@build-test.barrywardell.net"
msg['To'] = "carpetx-developers@einsteintoolkit.org"

# Record the MIME types of both parts - text/plain and text/html.
part1 = MIMEText(text, 'plain')
part2 = MIMEText(html, 'html')

# Attach parts into message container.
# According to RFC 2046, the last part of a multipart message, in this case
# the HTML message, is best and preferred.
msg.attach(part1)
msg.attach(part2)
 
# Send the message via our own SMTP server.
s = smtplib.SMTP('mail.einsteintoolkit.org')
s.send_message(msg)
s.quit()
