from gNB.crc import crc_encoder


data = '1101011111'
key = "10011"

data = [int(s) for s in data]
key = [int(s) for s in key]
print(data)

ans = crc_encoder(data, key)
print(ans)
