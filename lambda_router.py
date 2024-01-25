import subprocess as sb
import sys
import json as JSON

def run_crispresso(command):
    try:
        return_code = sb.call(command, shell=True)
        return return_code
    except sb.CalledProcessError as e:
        return e.output.decode(sys.stdout.encoding)


def lambda_handler(event, context):
    try:
        prepend_path = 'PATH=$PATH:/var/task/micromamba/bin '

        if 'body' in event.keys():
            event_dict = JSON.loads(event['body'])
            command = event_dict['command']
        else:
            command = event['command']

        commands = command.split(' && ')
        output = run_crispresso(prepend_path + commands[0])
        curl_command = str(commands[1])
        curl_command = curl_command.replace('"$?"', str(output))
        print(curl_command)
        output = run_crispresso(curl_command)

        # Create a response object with the output
        status = 200 if output == 0 else 500
        return JSON.dumps({
            "statusCode": status,
            "return_code": output,
            'command': command,
            'event': event
        })
    except Exception as e:
        return JSON.dumps({"statusCode": 401, "body": {"error": str(e)}, "return_code": -1, "command": command, "event": event})