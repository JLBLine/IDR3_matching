log_names = open('../IDR3multi_reject-eyeball_sources_v2.txt').read().split('\n')

log_names = [name for name in log_names if name!='']

old_choices = open('IDR3_matching_log_v1.txt').read().split('\n')

replace_lines = []

for i in xrange(len(old_choices)):
	line_entry = ''
	line = old_choices[i]
	if '#' in line:
		pass
	else:
		if line[:3] == 'GLM':
			line_entry += line
			line_entry += '\n'
			try:
				if 'GLM' not in old_choices[i+1]:
					line_entry += old_choices[i+1]
					line_entry += '\n'
					if 'GLM' not in old_choices[i+2]:
						line_entry += old_choices[i+2]
						line_entry += '\n'
						if 'GLM' not in old_choices[i+3]:
							line_entry += old_choices[i+3]
							line_entry += '\n'
							if 'GLM' not in old_choices[i+4]:
								line_entry += old_choices[i+4]
								line_entry += '\n'
								if 'GLM' not in old_choices[i+5]:
									line_entry += old_choices[i+5]
									line_entry += '\n'
									if 'GLM' not in old_choices[i+6]:
										line_entry += old_choices[i+6]
										line_entry += '\n'
			except:
				pass
	replace_lines.append(line_entry)
	
print len(old_choices),len(replace_lines)
	
new_log = open('IDR3_matching_log_v2.txt','w+')

if 'GLM000038+1213' in new_log:
	print 'up'
	
for name in log_names:
	written = 'no'
	
	for line in replace_lines:
		if name in line:
			if name == 'GLM000038+1213': print 'yup'
			#print name,line
			if 'GLM000038+1213' in line:
				print name,line
			
			new_log.write(line)
			written = 'yes'
			
	if written == 'no': new_log.write(name+'\n')
	
	
#for line in replace_lines:
	#new_log.write(line)
new_log.close()

sources = open('../IDR3multi_reject-eyeball_sources_v2.txt').read().split('\n')

choices = open('IDR3multi_reject-eyeball_choices_v1.txt').read().split('\n')

new_choices = open('IDR3multi_reject-eyeball_choices_v2.txt','w+')

for name in sources:
	written = 'no'
	for choice in choices:
		if name in choice:
			new_choices.write(choice+'\n')
			written = 'yes'
	if written == 'no':
		new_choices.write(name+'\n')
		
new_choices.close()