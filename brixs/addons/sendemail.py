#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Send email function"""

from getpass import getpass
import smtplib

def send_email(recipient, subject, body, sender, password=None):
    """Send email.

    If using a gmail sender email turn ``Allow less secure apps`` to ON.
    Be aware that this makes it easier for others to gain access to your account.

    Warning:
        Be aware that writing a password on a python script is a bad practice.
        It is recommended creating a development email just for usage with this
        function and never share sensitive information via this account.

    Args:
        recipient (str): email of the recipient.
        subject (str): email subject.
        body (str): email body.
        sender (str): sender email.
        password (str, optional): password of the sender email. If ``None``, it will
            securely ask the user to input the password.

    Returns:
        None
    """
    # assert smtplibok, 'send_email() cannot send email\nError: python package `smtplib` not found\nmaybe install it via ``pip install smtplib``' 
    if password is None:
        # assert getpassok, 'send_email() cannot ask for password securely\nError: python package `getpass` not found\nmaybe install it via ``pip install getpass``' 
        password = getpass()
    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.ehlo()
    server.login(sender, password)

    # Enter the headers of the email
    headers = '\r\n'.join(['from: ' + sender
                           , 'subject: ' + subject
                           , 'to: ' + recipient
                           , 'mime-version: 1.0'
                           , 'content-type: text/html'
                           ])
    # Enter the text of the body of the email
    body_of_email = body

    # Tie the headers and body together into the email's content
    content = headers + '\r\n\r\n' + body_of_email

    server.sendmail(sender, recipient, content)
    server.close()
    return

